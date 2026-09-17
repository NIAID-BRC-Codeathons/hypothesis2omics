"""Prototype: per-study GEO-readiness gate, run after ImmPort search and before
full retrieval (hypothesis2omics's "step 4", immport_fetch_module.fetch_immport_datasets).

For each candidate SDY accession, downloads only the small *_Tab.zip release file
(a few MB, not the full study) and inspects its contents:

    ready  = expsample_public_repository.txt is present AND has at least one row
             with REPOSITORY_NAME == "GEO" and a non-empty REPOSITORY_ACCESSION.

study_link.txt is recorded if present (its NAME/VALUE rows -- e.g. a
clinicaltrials.gov link) but never gates the verdict: it is a study-level table of
links to *any* external resource, not GEO-specific, and a study can have one without
having GEO-linked sample data (SDY1479), or have GEO-linked sample data without one
(SDY63).

This is a PROTOTYPE living outside hypothesis2omics on purpose. It imports the real
manifest/download functions from that repo's data/immport_fetch_module.py rather than
reimplementing them, so there is one source of truth for how ImmPort is queried; it
writes nothing back into that repo.

Run:
    python3 study_readiness_module.py SDY63 SDY1479 \
        --hypothesis2omics-root ~/Desktop/hypothesis2omics \
        --api-key-file ~/Desktop/h2o-scratch/.secrets/immport-key.json \
        --output-dir runs/READINESS01
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import sys
import time
import zipfile
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logger = logging.getLogger("study_readiness")

MODULE_VERSION = "0.1.0"
GEO_REPOSITORY_NAME = "geo"


@dataclass
class ReadinessRecord:
    sdy_id: str
    tab_zip_path: str | None
    has_expsample_public_repository: bool
    geo_linked_sample_count: int
    has_study_link: bool
    study_link_targets: list[str] = field(default_factory=list)
    ready: bool = False
    reason: str = ""

    def to_row(self) -> dict[str, Any]:
        row = asdict(self)
        row["study_link_targets"] = "; ".join(self.study_link_targets)
        return row


ROW_FIELDS = [
    "sdy_id",
    "tab_zip_path",
    "has_expsample_public_repository",
    "geo_linked_sample_count",
    "has_study_link",
    "study_link_targets",
    "ready",
    "reason",
]


def _load_immport_fetch_module(hypothesis2omics_root: Path):
    """Import the real fetch module from the repo by path, without editing it."""
    data_dir = str((hypothesis2omics_root / "data").resolve())
    if data_dir not in sys.path:
        sys.path.insert(0, data_dir)
    import immport_fetch_module as m  # noqa: PLC0415

    return m


def _parse_expsample_public_repository(content: bytes) -> tuple[bool, int]:
    reader = csv.DictReader(io.StringIO(content.decode("utf-8", errors="replace")), delimiter="\t")
    fieldnames = [f.strip().upper() for f in (reader.fieldnames or [])]
    if "REPOSITORY_NAME" not in fieldnames or "REPOSITORY_ACCESSION" not in fieldnames:
        return True, 0
    count = 0
    for row in reader:
        name = str(row.get("REPOSITORY_NAME") or "").strip().lower()
        accession = str(row.get("REPOSITORY_ACCESSION") or "").strip()
        if name == GEO_REPOSITORY_NAME and accession:
            count += 1
    return True, count


def _parse_study_link(content: bytes) -> list[str]:
    reader = csv.DictReader(io.StringIO(content.decode("utf-8", errors="replace")), delimiter="\t")
    targets = []
    for row in reader:
        name = str(row.get("NAME") or row.get("name") or "").strip()
        value = str(row.get("VALUE") or row.get("value") or "").strip()
        if name or value:
            targets.append(f"{name} -> {value}" if name else value)
    return targets


def _find_member(namelist: list[str], filename: str) -> str | None:
    matches = [n for n in namelist if n.rsplit("/", 1)[-1] == filename]
    return matches[0] if matches else None


def check_study_readiness(sdy_id: str, session, cache_dir: Path, fetch_module) -> ReadinessRecord:
    entries = fetch_module._fetch_filepath_manifest(sdy_id, session)
    tab_entry = next((e for e in entries if str(e.get("path", "")).endswith("_Tab.zip")), None)
    if tab_entry is None:
        return ReadinessRecord(
            sdy_id=sdy_id,
            tab_zip_path=None,
            has_expsample_public_repository=False,
            geo_linked_sample_count=0,
            has_study_link=False,
            ready=False,
            reason="no *_Tab.zip release file in ImmPort manifest",
        )

    sdy_dir = cache_dir / sdy_id
    relative_parts = Path(tab_entry["path"].lstrip("/")).parts
    if relative_parts and relative_parts[0].upper() == sdy_id.upper():
        relative_parts = relative_parts[1:]
    destination = sdy_dir.joinpath(*relative_parts)

    artifact = fetch_module._download_file(
        tab_entry["path"],
        tab_entry["fileUUID"],
        destination,
        session,
        force=False,
        artifact_type="release_file",
        expected_size=tab_entry.get("filesizeBytes"),
    )
    if artifact.status not in ("success", "cached"):
        return ReadinessRecord(
            sdy_id=sdy_id,
            tab_zip_path=str(destination),
            has_expsample_public_repository=False,
            geo_linked_sample_count=0,
            has_study_link=False,
            ready=False,
            reason=f"Tab.zip download failed: {artifact.error}",
        )

    with zipfile.ZipFile(destination) as z:
        namelist = z.namelist()
        expsample_member = _find_member(namelist, "expsample_public_repository.txt")
        study_link_member = _find_member(namelist, "study_link.txt")

        has_expsample = False
        geo_count = 0
        if expsample_member:
            with z.open(expsample_member) as f:
                has_expsample, geo_count = _parse_expsample_public_repository(f.read())

        has_study_link = study_link_member is not None
        study_link_targets: list[str] = []
        if study_link_member:
            with z.open(study_link_member) as f:
                study_link_targets = _parse_study_link(f.read())

    if not has_expsample:
        reason = "no expsample_public_repository.txt in Tab.zip"
    elif geo_count == 0:
        reason = "expsample_public_repository.txt has no GEO rows"
    else:
        reason = f"{geo_count} GEO-linked sample(s) in expsample_public_repository.txt"

    return ReadinessRecord(
        sdy_id=sdy_id,
        tab_zip_path=str(destination),
        has_expsample_public_repository=has_expsample,
        geo_linked_sample_count=geo_count,
        has_study_link=has_study_link,
        study_link_targets=study_link_targets,
        ready=has_expsample and geo_count > 0,
        reason=reason,
    )


def assess_studies(
    sdy_ids: list[str],
    hypothesis2omics_root: Path,
    output_dir: Path,
    api_key_file: str | None = None,
    cache_dir: Path | None = None,
) -> dict[str, Any]:
    started = time.monotonic()
    fetch_module = _load_immport_fetch_module(hypothesis2omics_root)
    session = fetch_module.ImmportSession(api_key_file=api_key_file)
    cache = cache_dir or (output_dir / "immport_cache")
    cache.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = [check_study_readiness(sdy, session, cache, fetch_module) for sdy in sdy_ids]

    tsv_path = output_dir / "study_readiness.tsv"
    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_FIELDS, delimiter="\t")
        writer.writeheader()
        for r in records:
            writer.writerow(r.to_row())

    ready_ids = [r.sdy_id for r in records if r.ready]
    excluded = [r for r in records if not r.ready]

    (output_dir / "ready_studies.txt").write_text("\n".join(ready_ids) + ("\n" if ready_ids else ""))

    excluded_path = output_dir / "excluded_studies.tsv"
    with excluded_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_FIELDS, delimiter="\t")
        writer.writeheader()
        for r in excluded:
            writer.writerow(r.to_row())

    provenance = {
        "run_started_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool": "study_readiness_module",
        "tool_version": MODULE_VERSION,
        "studies_checked": len(records),
        "ready": ready_ids,
        "excluded": [r.sdy_id for r in excluded],
        "outputs": {
            "study_readiness_tsv": str(tsv_path),
            "ready_studies_txt": str(output_dir / "ready_studies.txt"),
            "excluded_studies_tsv": str(excluded_path),
        },
        "duration_sec": round(time.monotonic() - started, 3),
    }
    (output_dir / "study_readiness.provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    return {"records": records, "provenance": provenance}


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Gate ImmPort study accessions on GEO-linked sample data before retrieval."
    )
    parser.add_argument("sdy_ids", nargs="+", help="ImmPort study accessions to check")
    parser.add_argument("--hypothesis2omics-root", type=Path, required=True,
                        help="Path to a checkout of the hypothesis2omics repo (read-only import)")
    parser.add_argument("--api-key-file", required=True, help="ImmPort API key JSON file")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, default=None,
                        help="Where to cache downloaded Tab.zip files (default: <output-dir>/immport_cache)")
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    result = assess_studies(
        args.sdy_ids,
        hypothesis2omics_root=args.hypothesis2omics_root,
        output_dir=args.output_dir,
        api_key_file=args.api_key_file,
        cache_dir=args.cache_dir,
    )
    for r in result["records"]:
        logger.info("%s ready=%s (%s)", r.sdy_id, r.ready, r.reason)
    logger.info("%s/%s studies ready", len(result["provenance"]["ready"]), len(result["records"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
