# /// script
# requires-python = ">=3.11"
# dependencies = ["mcp>=2.0", "httpx>=0.27", "pyyaml>=6"]
# ///
"""Mini MCP server: keywords / a search spec -> ImmPort study accession numbers.

Wraps exactly one ImmPort endpoint: GET /api/search/study (study_search).
No auth required for this endpoint.
"""

import sys
from typing import Any

import httpx
import yaml
from mcp.server.mcpserver import MCPServer

API = "https://www.immport.org/data/query/api/search/study"
FIELDS = "study_accession,brief_title,condition_studied,research_focus"

# Facet query params of study_search. Values are controlled vocabulary except the
# accession/age ones; use list_facet_values() to get the valid values.
FACETS = frozenset(
    "studyAccession subjectAccession programName conditionOrDisease researchFocus "
    "clinicalTrial sex race ethnicity species minAge ageRange assayMethod "
    "biosampleType hasAssessment hasLabTest searchFields".split()
)

server = MCPServer("immport-search")


def _get(client: httpx.Client, params: dict) -> dict:
    r = client.get(API, params=params, timeout=60)
    r.raise_for_status()
    return r.json()


def _filters(raw: dict[str, Any]) -> dict[str, str]:
    bad = set(raw) - FACETS
    if bad:
        raise ValueError(f"unknown facet(s) {sorted(bad)}; valid: {sorted(FACETS)}")
    # ImmPort ORs values within one facet, ANDs across facets.
    return {k: ",".join(map(str, v)) if isinstance(v, (list, tuple)) else str(v) for k, v in raw.items()}


def _hits(client: httpx.Client, term: str | None, filters: dict, page_size: int, fields: str) -> list[dict]:
    params = {"pageSize": page_size, "sourceFields": fields, **filters}
    if term:
        params["term"] = term
    return _get(client, params)["hits"]["hits"]


def _search(keywords: list[str], limit: int, exact: bool, filters: dict | None = None) -> dict:
    if not keywords:
        raise ValueError("keywords must not be empty")
    limit = max(1, min(int(limit), 1000))  # API caps pageSize at 1000
    # term ANDs its words, so one query per keyword; results are unioned.
    term = " ".join(f'"{k}"' if exact else k for k in keywords)
    with httpx.Client() as c:
        params = {"term": term, "pageSize": limit, "sourceFields": FIELDS, **_filters(filters or {})}
        d = _get(c, params)
    return {
        "total": d["hits"]["total"]["value"],
        "studies": [
            {
                "accession": h["_id"],
                "score": h.get("_score"),
                **{k: v for k, v in h.get("_source", {}).items() if k != "study_accession"},
            }
            for h in d["hits"]["hits"]
        ],
    }


def _spec_search(spec: dict) -> list[str]:
    groups: dict[str, list[str]] = spec.get("search_terms") or {}
    if not isinstance(groups, dict):
        raise ValueError("search_terms must be a mapping of group name -> list of terms")
    filters = _filters(spec.get("filters") or {})
    require = set(spec.get("require") or [])
    unknown = require - set(groups)
    if unknown:
        raise ValueError(f"require names unknown group(s) {sorted(unknown)}")
    min_groups = int(spec.get("min_groups", len(require) or 1))
    limit = max(1, min(int(spec.get("limit", 100)), 1000))
    per_term = max(1, min(int(spec.get("per_term", 200)), 1000))

    matched: dict[str, set[str]] = {}  # accession -> group names
    scores: dict[str, float] = {}
    with httpx.Client() as c:
        if not groups:  # facet-only search
            for h in _hits(c, None, filters, limit, "study_accession"):
                matched[h["_id"]], scores[h["_id"]] = {"*"}, h.get("_score") or 0.0
            min_groups = 1
        for group, terms in groups.items():
            for term in dict.fromkeys(terms or []):  # dedupe, keep order
                for h in _hits(c, str(term), filters, per_term, "study_accession"):
                    matched.setdefault(h["_id"], set()).add(group)
                    scores[h["_id"]] = scores.get(h["_id"], 0.0) + (h.get("_score") or 0.0)

    hits = [a for a, g in matched.items() if len(g) >= min_groups and require <= g]
    hits.sort(key=lambda a: (-len(matched[a]), -scores[a], a))
    return hits[:limit]


@server.tool()
def search_studies(
    keywords: list[str], limit: int = 10, exact_phrases: bool = False, filters: dict[str, Any] | None = None
) -> dict:
    """Convert keywords to ImmPort study accession numbers (SDY...).

    Args:
        keywords: search terms. NOTE: ImmPort ANDs all words of a query, so
            ["influenza", "malaria"] means influenza AND malaria (usually 0 hits).
            For OR / multi-facet logic use search_spec instead.
        limit: max studies to return (1-1000).
        exact_phrases: quote each keyword so only the exact phrase matches.
        filters: optional facet filters, e.g. {"species": "Homo sapiens",
            "assayMethod": ["RNA sequencing", "Transcription profiling by array"]}.

    Returns: {"total": <total matches>, "studies": [{"accession", "score", ...}]}
    """
    return _search(keywords, limit, exact_phrases, filters)


@server.tool()
def search_spec(spec: str | dict[str, Any]) -> list[str]:
    """Run a multi-facet keyword spec against ImmPort; returns accession IDs only.

    Each term is queried separately (ImmPort has no OR operator), so terms within a
    group are ORed and groups are scored by how many of them a study matches.
    Results are sorted best-first: most groups matched, then summed relevance.

    Args:
        spec: YAML or JSON (string or object):

            repositories: [ImmPort]        # optional, ignored (ImmPort only)
            search_terms:                  # group name -> terms (OR within a group)
              intervention: [YF-17D, yellow fever vaccine]
              assay: [RNA-seq, transcriptomics]
            require: [intervention]        # optional: groups a study MUST match
            min_groups: 2                  # optional: min groups matched (default:
                                           # len(require) or 1)
            filters:                       # optional facet filters, ANDed with every
              species: Homo sapiens        # term query. Controlled vocabulary -
              assayMethod: [RNA sequencing]  # see list_facet_values()
            limit: 100                     # optional, max accessions returned
            per_term: 200                  # optional, hits fetched per term

        Free-text keys such as required_features are ignored - express hard
        requirements as `require` groups or `filters`.

    Returns: list of study accessions (e.g. ["SDY1264", "SDY1289"]), best first.
    """
    parsed = yaml.safe_load(spec) if isinstance(spec, str) else spec
    if not isinstance(parsed, dict):
        raise ValueError("spec must be a YAML/JSON mapping")
    return _spec_search(parsed)


@server.tool()
def list_facet_values(facet: str | None = None) -> dict[str, list[str]]:
    """List the controlled-vocabulary values accepted by `filters`.

    Args:
        facet: optional facet name (e.g. "assayMethod") to return just that one.

    Returns: {facet name: [values, ...]} with the most-used values first.
    """
    with httpx.Client() as c:
        aggs = _get(c, {"pageSize": 0})["aggregations"]
    # Aggregation keys are snake_case; facet params are camelCase.
    out = {
        "".join(w.title() if i else w for i, w in enumerate(k.split("_"))): [
            b["key"] for b in v["buckets"] if isinstance(b, dict)
        ]
        for k, v in aggs.items()
        if isinstance(v, dict) and isinstance(v.get("buckets"), list)
    }
    out["sex"] = out.pop("gender", [])
    out = {k: v for k, v in out.items() if k in FACETS}
    return {facet: out[facet]} if facet else out


def check() -> None:
    out = _search(["malaria"], 3, False)
    assert out["total"] > 0, out
    assert len(out["studies"]) == 3, out
    assert all(s["accession"].startswith("SDY") for s in out["studies"]), out

    vocab = list_facet_values()
    assert "Homo sapiens" in vocab["species"], vocab["species"][:5]
    assert "RNA sequencing" in vocab["assayMethod"], vocab["assayMethod"][:5]

    spec = """
    search_terms:
      intervention: [YF-17D, YF17D, yellow fever vaccine, yellow fever vaccination]
      predictor: [GCN2, EIF2AK4]
      outcome: [CD8, CD8+ T cell, T-cell response, vaccine response]
      assay: [transcriptomics, RNA-seq, microarray, gene expression]
    require: [intervention]
    filters:
      species: Homo sapiens
    limit: 10
    """
    accs = search_spec(spec)
    assert accs and all(a.startswith("SDY") for a in accs), accs
    narrow = search_spec({"search_terms": {"i": ["YF-17D"]}, "min_groups": 1, "limit": 5})
    assert set(narrow) <= set(_spec_search({"search_terms": {"i": ["YF-17D"]}, "limit": 1000})), narrow
    print("ok", out["studies"][0]["accession"], accs)


if __name__ == "__main__":
    check() if "--check" in sys.argv else server.run()
