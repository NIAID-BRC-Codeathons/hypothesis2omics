"""Step 4: per-study GEO-readiness gate, run after ImmPort search (keyword-to-immport's
search_spec) and before full retrieval (data/immport_fetch_module.fetch_immport_datasets).

For each candidate SDY accession, downloads only the small *_Tab.zip release file
(a few MB, not the full study) and runs the real ImmPort sample-link parser
(data/immport_parse_module.py) against it -- the same parser step 6
(immport_batch_parse.py) would otherwise run separately:

    ready  = the 8-table ImmPort join succeeds AND yields at least one sample row
             with repository_name == "GEO" and a non-empty repository_accession.

Because this reuses the real parser rather than checking for one required file by
name, "ready" now also catches any of the other 7 required tables being missing or
malformed -- failure modes step 6 would otherwise hit later, on a study step 4 had
already waved through. One parse, two outputs: the readiness verdict, and every ready
study's contribution to a combined sample_manifest.tsv (same shape
immport_batch_parse.py's combined output already uses), written alongside the other
readiness artifacts. Not-ready studies never contribute rows to it.

study_link.txt is recorded if present (its NAME/VALUE rows -- e.g. a
clinicaltrials.gov link) but never gates the verdict: it is a study-level table of
links to *any* external resource, not GEO-specific, and a study can have one without
having GEO-linked sample data (SDY1479), or have GEO-linked sample data without one
(SDY63).

Imports the real manifest/download and parse functions from
../../data/immport_fetch_module.py and ../../data/immport_parse_module.py by path,
rather than reimplementing them, so there is one source of truth for how ImmPort is
queried and joined. The parse side needs pandas -- run this under an environment that
has it (e.g. `uv run` against the repo root's pyproject.toml), not a bare venv.

Run:
    python3 study_readiness_module.py SDY63 SDY1479 \
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

MODULE_VERSION = "0.2.0"
GEO_REPOSITORY_NAME = "geo"


@dataclass
class ReadinessRecord:
    sdy_id: str
    tab_zip_path: str | None
    parsed_successfully: bool
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
    "parsed_successfully",
    "geo_linked_sample_count",
    "has_study_link",
    "study_link_targets",
    "ready",
    "reason",
]


def _data_dir() -> str:
    return str(Path(__file__).resolve().parents[2] / "data")


def _load_immport_fetch_module():
    """Import the real fetch module by path, the way build_test_spec.py imports
    keyword-to-immport/server.py -- one copy, resolved from this file's location
    rather than asked for on the command line."""
    data_dir = _data_dir()
    if data_dir not in sys.path:
        sys.path.insert(0, data_dir)
    import immport_fetch_module as m  # noqa: PLC0415

    return m


def _load_immport_parse_module():
    """Import the real sample-link parser by path, same pattern as
    _load_immport_fetch_module. This is what turns the readiness check into the
    authoritative parse: one join, reused for both the verdict and the manifest,
    instead of a lighter hand-rolled check here plus a full reparse in step 6."""
    data_dir = _data_dir()
    if data_dir not in sys.path:
        sys.path.insert(0, data_dir)
    import immport_parse_module as m  # noqa: PLC0415

    return m


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


def check_study_readiness(
    sdy_id: str, session, cache_dir: Path, fetch_module, parse_module
) -> tuple[ReadinessRecord, Any | None]:
    """Returns the readiness verdict, and -- only when ready -- the study's rows
    for the combined sample_manifest.tsv (a pandas DataFrame, or None)."""
    entries = fetch_module._fetch_filepath_manifest(sdy_id, session)
    tab_entry = next((e for e in entries if str(e.get("path", "")).endswith("_Tab.zip")), None)
    if tab_entry is None:
        return ReadinessRecord(
            sdy_id=sdy_id,
            tab_zip_path=None,
            parsed_successfully=False,
            geo_linked_sample_count=0,
            has_study_link=False,
            ready=False,
            reason="no *_Tab.zip release file in ImmPort manifest",
        ), None

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
            parsed_successfully=False,
            geo_linked_sample_count=0,
            has_study_link=False,
            ready=False,
            reason=f"Tab.zip download failed: {artifact.error}",
        ), None

    # study_link.txt is informational only, unrelated to the 8-table parse below --
    # a separate, lightweight pass over the same local zip.
    with zipfile.ZipFile(destination) as z:
        study_link_member = _find_member(z.namelist(), "study_link.txt")
        has_study_link = study_link_member is not None
        study_link_targets: list[str] = []
        if study_link_member:
            with z.open(study_link_member) as f:
                study_link_targets = _parse_study_link(f.read())

    try:
        frame = parse_module.parse_immport_sample_links(destination)
    except parse_module.ImmportParseError as exc:
        return ReadinessRecord(
            sdy_id=sdy_id,
            tab_zip_path=str(destination),
            parsed_successfully=False,
            geo_linked_sample_count=0,
            has_study_link=has_study_link,
            study_link_targets=study_link_targets,
            ready=False,
            reason=f"ImmPort sample-link parse failed: {exc}",
        ), None

    geo_rows = frame[frame["repository_name"].str.lower() == GEO_REPOSITORY_NAME]
    geo_count = int((geo_rows["repository_accession"].str.strip() != "").sum())

    if geo_count == 0:
        reason = "sample manifest parsed but has no GEO-linked rows"
    else:
        reason = f"{geo_count} GEO-linked sample(s) in the parsed sample manifest"

    record = ReadinessRecord(
        sdy_id=sdy_id,
        tab_zip_path=str(destination),
        parsed_successfully=True,
        geo_linked_sample_count=geo_count,
        has_study_link=has_study_link,
        study_link_targets=study_link_targets,
        ready=geo_count > 0,
        reason=reason,
    )
    return record, (frame if geo_count > 0 else None)


def assess_studies(
    sdy_ids: list[str],
    output_dir: Path,
    api_key_file: str | None = None,
    cache_dir: Path | None = None,
) -> dict[str, Any]:
    started = time.monotonic()
    fetch_module = _load_immport_fetch_module()
    parse_module = _load_immport_parse_module()
    session = fetch_module.ImmportSession(api_key_file=api_key_file)
    cache = cache_dir or (output_dir / "immport_cache")
    cache.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    records: list[ReadinessRecord] = []
    manifest_frames = []
    for sdy in sdy_ids:
        record, frame = check_study_readiness(sdy, session, cache, fetch_module, parse_module)
        records.append(record)
        if frame is not None:
            manifest_frames.append(frame)

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

    manifest_path = None
    manifest_rows = 0
    if manifest_frames:
        import pandas as pd  # noqa: PLC0415 -- only needed on this path

        combined = pd.concat(manifest_frames, ignore_index=True)
        combined = combined.loc[:, parse_module.OUTPUT_COLUMNS].sort_values(
            ["study_accession", "expsample_accession"]
        )
        manifest_path = output_dir / "sample_manifest.tsv"
        combined.to_csv(manifest_path, sep="\t", index=False, lineterminator="\n")
        manifest_rows = len(combined)

    provenance = {
        "run_started_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool": "study_readiness_module",
        "tool_version": MODULE_VERSION,
        "studies_checked": len(records),
        "ready": ready_ids,
        "excluded": [r.sdy_id for r in excluded],
        "sample_manifest": {
            "path": str(manifest_path) if manifest_path else None,
            "rows": manifest_rows,
            "studies_included": [r.sdy_id for r in records if r.ready],
        },
        "outputs": {
            "study_readiness_tsv": str(tsv_path),
            "ready_studies_txt": str(output_dir / "ready_studies.txt"),
            "excluded_studies_tsv": str(excluded_path),
            "sample_manifest_tsv": str(manifest_path) if manifest_path else None,
        },
        "duration_sec": round(time.monotonic() - started, 3),
    }
    (output_dir / "study_readiness.provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    return {"records": records, "provenance": provenance}


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Gate ImmPort study accessions on GEO-linked sample data before "
                    "retrieval, and produce a combined sample_manifest.tsv for the "
                    "studies that pass."
    )
    parser.add_argument("sdy_ids", nargs="+", help="ImmPort study accessions to check")
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
        output_dir=args.output_dir,
        api_key_file=args.api_key_file,
        cache_dir=args.cache_dir,
    )
    for r in result["records"]:
        logger.info("%s ready=%s (%s)", r.sdy_id, r.ready, r.reason)
    logger.info("%s/%s studies ready", len(result["provenance"]["ready"]), len(result["records"]))
    manifest = result["provenance"]["sample_manifest"]
    if manifest["path"]:
        logger.info("wrote %s rows to %s", manifest["rows"], manifest["path"])
    else:
        logger.info("no ready studies, no sample_manifest.tsv written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
