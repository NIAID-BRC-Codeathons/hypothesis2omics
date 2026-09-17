#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["mcp>=2.0", "openai>=1.40", "httpx>=0.27", "pyyaml>=6"]
# ///
"""MCP server: a plain-text hypothesis -> components -> test spec -> ImmPort search spec.

Steps 1-3 of the pipeline, exposed as tools. The CLIs remain the reference
implementation; this imports them rather than reimplementing them, so there is one
copy of every rule.

    parse_hypothesis(text, outdir)          STEP 1  (Argo)  -> 01_parsed.yaml
            |
            |   <- a person reads 01_parsed.yaml here and edits what is wrong
            v
    build_test_spec(parsed_path)            STEP 2  (Argo)  -> 02_test_spec.json
                                            STEP 3  (local) -> 03_search_spec.yaml
    derive_search_spec(test_spec_path)      STEP 3 again, after a hand edit
    audit_search_terms(search_spec_path)    per-term hit counts, ImmPort reads only

ON THE HUMAN GATE. build_test_spec takes a filesystem path, never an inline object,
so 01_parsed.yaml must exist before step 2 can run. That is a speed bump, not an
enforcement: an agent can write the file and call step 2 in the same turn. Provenance
records the input file's sha256 and makes no claim that anybody read it. If you need
evidence of review, it has to come from outside this server.

Every tool writes its numbered artifact to disk and also returns the content, so an
MCP run leaves the same audit trail a CLI run does.

Errors are raised, not swallowed. An empty Argo completion, a hand edit that breaks
the parsed schema, or a spec the validator rejects surfaces as a failed tool call --
never as a plausible-looking empty result.

Notes to stderr (model fallback, budget growth, repair passes) are safe: MCP speaks
over stdout, and nothing in this pipeline writes there outside of --check.

Install:
    pip install mcp openai httpx pyyaml

Run:
    ARGO_USER=ac.yourname python3 hypothesis_mcp.py --check   # smoke test, no client
    ARGO_USER=ac.yourname python3 hypothesis_mcp.py           # serve over stdio

Requires a connection to the Argonne-auth network.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

import yaml
from mcp.server.mcpserver import MCPServer

# MCP clients launch a server from an arbitrary working directory, so the sibling
# modules have to be findable by path rather than by cwd.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import parse_hypothesis as step1  # noqa: E402
from h2o_common import ARGO_BASE_URL, ARGO_USER, DEFAULT_MODEL  # noqa: E402

server = MCPServer("hypothesis2omics")

# Per-attempt ceiling on an Argo call. Step 2 is normally 30-90s; the OpenAI client
# retries on its own, so a hard hang costs a small multiple of this.
ARGO_TIMEOUT_S = float(os.environ.get("H2O_TIMEOUT_S", "100"))

DEFAULT_MAX_TOKENS = 4000

# What build_test_spec.py writes into the step 3 header. Repeated here so
# derive_search_spec produces a byte-identical file when it re-derives one.
STEP3_HEADER = {
    "step": "3-search-spec",
    "derived_from": "02_test_spec.json",
    "method": "deterministic collection of components[*].search_terms",
}


def _argo() -> Any:
    """An Argo client with a timeout. h2o_common.client() does not take one."""
    from openai import OpenAI

    return OpenAI(api_key=ARGO_USER, base_url=ARGO_BASE_URL, timeout=ARGO_TIMEOUT_S)


def _step2() -> Any:
    """Import steps 2-3 lazily.

    build_test_spec.py pulls in server.py, httpx and the ImmPort layer. Deferring it
    to first use means a client that only ever calls parse_hypothesis keeps step 1's
    property of running with the repository unreachable.
    """
    import build_test_spec

    return build_test_spec


def _outdir(path: str) -> Path:
    out = Path(path).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    return out


def _read(path: str, what: str) -> Path:
    src = Path(path).expanduser()
    if not src.is_file():
        raise FileNotFoundError(f"no {what} at {src}")
    return src


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


@server.tool()
def parse_hypothesis(
    hypothesis_text: str, outdir: str, model: str = DEFAULT_MODEL,
    max_tokens: int = DEFAULT_MAX_TOKENS
) -> dict:
    """Step 1: decompose a plain-English hypothesis into its named components.

    Writes 01_parsed.yaml, which is the entire interface to step 2. One Argo call,
    typically 10-30s. Touches neither ImmPort nor the network beyond Argo.

    Args:
        hypothesis_text: the hypothesis, in prose. One sentence is normal.
        outdir: directory to write 01_parsed.yaml into; created if absent.
        model: Argo model ID (not the vendor's). Falls back to claudesonnet5 if the
            requested model returns an empty completion; provenance records which
            model actually answered.
        max_tokens: text budget for the call.

    Returns: {"path", "parsed", "warnings", "yaml", "provenance", "next_step"}.

        `warnings` is the parser naming what it inferred rather than read -- an
        inferred cell_type, an outcome with no measurable quantity. Nothing
        downstream branches on it. Surface it to whoever is reviewing the file.

    Raises: EmptyCompletion if Argo returns no text on either model; ValueError if
        the model cannot produce a parse the validator accepts.
    """
    out = _outdir(outdir)
    parsed, provenance = step1.parse_hypothesis(_argo(), hypothesis_text, model, max_tokens)
    text = step1.parsed_to_yaml(parsed, hypothesis_text, provenance)
    path = out / "01_parsed.yaml"
    path.write_text(text)
    return {
        "path": str(path),
        "parsed": {f: parsed.get(f) for f in step1.PARSED_FIELDS},
        "warnings": parsed.get("warnings") or [],
        "yaml": text,
        "provenance": provenance,
        "next_step": f"review {path.name}, correct anything wrong, then call "
                     f"build_test_spec with parsed_path={str(path)!r}",
    }


@server.tool()
def build_test_spec(
    parsed_path: str, outdir: str | None = None, model: str = DEFAULT_MODEL,
    hypothesis_id: str | None = None, max_tokens: int = DEFAULT_MAX_TOKENS
) -> dict:
    """Steps 2-3: an edited 01_parsed.yaml -> a test spec and an ImmPort search spec.

    Step 2 expands each component into synonyms, candidate measurements and checkable
    eligibility criteria (one Argo call, 30-90s, plus a vocabulary read from ImmPort).
    Step 3 then collects the per-component search_terms deterministically -- terms are
    lifted from step 2, never invented here.

    Takes a path, not an object: 01_parsed.yaml has to exist on disk. Nothing in this
    tool verifies that a person read it.

    Args:
        parsed_path: path to 01_parsed.yaml, hand-edits and all.
        outdir: where to write 02 and 03. Defaults to the parsed file's directory.
        model: Argo model ID.
        hypothesis_id: override the id the model generates.
        max_tokens: text budget for the call.

    Returns: {"test_spec_path", "search_spec_path", "test_spec", "search_spec",
        "search_spec_yaml", "provenance"}.

    Raises: ValueError listing every problem if the parsed file breaks its schema --
        a bad hand edit stops here rather than surfacing downstream; EmptyCompletion
        if Argo returns no text; ValueError if the spec fails validation after a
        repair pass.
    """
    src = _read(parsed_path, "parsed hypothesis (run parse_hypothesis first)")
    step2 = _step2()
    parsed, hypothesis = step1.load_parsed(src)

    spec, provenance = step2.build_test_spec(_argo(), hypothesis, parsed, model, max_tokens)
    if hypothesis_id:
        spec["hypothesis_id"] = hypothesis_id
    # Which bytes step 2 actually read. Not a claim that anyone reviewed them.
    provenance["parsed_file"] = src.name
    provenance["parsed_sha256"] = _digest(src)
    spec["provenance"] = provenance

    out = _outdir(outdir) if outdir else src.parent
    spec_path = out / "02_test_spec.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n")

    search = step2.to_search_spec(spec)
    search_yaml = step2.search_spec_to_yaml(search, False, STEP3_HEADER)
    search_path = out / "03_search_spec.yaml"
    search_path.write_text(search_yaml)

    return {
        "test_spec_path": str(spec_path),
        "search_spec_path": str(search_path),
        "test_spec": spec,
        "search_spec": search,
        "search_spec_yaml": search_yaml,
        "provenance": provenance,
    }


@server.tool()
def derive_search_spec(test_spec_path: str, poc_format: bool = False) -> dict:
    """Step 3 on its own: re-derive 03_search_spec.yaml from an edited 02_test_spec.json.

    Deterministic and offline apart from reading the file -- no Argo call, no ImmPort
    query. Use it after hand-editing search_terms, require or filters in the test spec,
    so the two files cannot drift apart.

    Args:
        test_spec_path: path to 02_test_spec.json.
        poc_format: emit only the Example_input_output.md keys, omitting the
            require/min_groups/filters block that actually constrains the query. A
            document rendering; not usable as input to a search.

    Returns: {"path", "search_spec", "yaml"}.
    """
    src = _read(test_spec_path, "test spec")
    step2 = _step2()
    spec = json.loads(src.read_text())
    search = step2.to_search_spec(spec)
    header = dict(STEP3_HEADER, derived_from=src.name)
    text = step2.search_spec_to_yaml(search, poc_format, header)
    path = src.parent / "03_search_spec.yaml"
    path.write_text(text)
    return {"path": str(path), "search_spec": search, "yaml": text}


@server.tool()
def audit_search_terms(search_spec_path: str) -> dict:
    """Report what each search term returns on its own, before trusting the union.

    One ImmPort query per term (no Argo call), so a spec with forty terms takes a
    while. `hits` is what a term returns alone; `unique` is what only that term found.

    A term with hits=0 contributes nothing. A term in a required group with a high
    `unique` is where false positives enter -- measured 2026-09-16, "17D" alone
    returned 14 studies of which 9 were unrelated, because short fragments of a
    compound identifier collide with arm and protocol codes.

    Args:
        search_spec_path: path to 03_search_spec.yaml.

    Returns: {"terms": [{"group", "term", "required_group", "hits", "unique"}, ...],
        "no_hits": [...], "widens_required_set": [...]}, best-first.
    """
    src = _read(search_spec_path, "search spec")
    step2 = _step2()
    spec = yaml.safe_load(src.read_text())
    if not isinstance(spec, dict):
        raise ValueError(f"{src} is not a YAML mapping")
    rows = step2.audit_terms(spec)
    return {
        "terms": rows,
        "no_hits": [r["term"] for r in rows if r["hits"] == 0],
        "widens_required_set": [r["term"] for r in rows
                                if r["required_group"] and r["unique"] >= 3],
    }


def check() -> None:
    hypothesis = ("GCN2/EIF2AK4 activity is associated with the magnitude of the CD8+ "
                  "T-cell response following YF-17D vaccination.")
    with tempfile.TemporaryDirectory() as tmp:
        one = parse_hypothesis(hypothesis, tmp)
        assert Path(one["path"]).is_file(), one
        assert one["parsed"]["outcome"], one["parsed"]
        assert one["parsed"]["predictor"], one["parsed"]
        assert isinstance(one["warnings"], list), one

        two = build_test_spec(one["path"], hypothesis_id="CHECK01")
        spec = two["test_spec"]
        assert spec["hypothesis_id"] == "CHECK01", spec["hypothesis_id"]
        assert spec["eligibility"]["scientific_criteria"], spec["eligibility"]
        assert len(spec["inclusion_criteria"]) == len(spec["eligibility"]["scientific_criteria"])
        assert spec["provenance"]["parsed_sha256"] == _digest(Path(one["path"]))
        assert "intervention" in two["search_spec"]["search_terms"], two["search_spec"]

        three = derive_search_spec(two["test_spec_path"])
        assert three["search_spec"] == two["search_spec"], "step 3 is not deterministic"

        # Audit a two-term spec rather than the real one: same code path, bounded cost.
        small = Path(tmp) / "small_search_spec.yaml"
        small.write_text("search_terms:\n  intervention:\n    - YF-17D\n    - nonsense-term\n"
                         "require:\n  - intervention\nlimit: 10\n")
        audit = audit_search_terms(str(small))
        assert len(audit["terms"]) == 2, audit["terms"]
        assert "nonsense-term" in audit["no_hits"], audit

    print("ok", one["provenance"]["model"], two["test_spec_path"],
          sorted(two["search_spec"]["search_terms"]))


if __name__ == "__main__":
    check() if "--check" in sys.argv else server.run()
