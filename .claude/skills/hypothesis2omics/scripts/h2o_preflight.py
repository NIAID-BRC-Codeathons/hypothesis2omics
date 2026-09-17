# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Preflight for the hypothesis2omics pipeline: can it run, where is it, was it run.

Purpose
    Answer three questions before an agent or a person calls a pipeline tool:
    (1) is the environment able to run the stage, (2) which stage is next, and
    (3) was every artifact on disk produced by the pipeline rather than by hand.

Inputs
    --repo PATH   checkout of hypothesis2omics (default: current directory)
    --json        emit the report as JSON instead of text

Outputs
    A report on stdout. Exit 0 = proceed, 1 = something named is blocking,
    2 = usage error.

External dependencies
    None. Standard library only, so it runs before `uv sync` has been done.
    Reads environment variables and the Loom/pi MCP config; writes nothing.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Literal

Level = Literal["block", "note", "ok"]

# Tools whose outputs are expected to carry a sibling *.provenance.json, versus
# stages that are known not to emit one. The second list is a documented gap,
# not a licence -- see references/file-contract.md.
PROVENANCE_EXEMPT_NOTE = "stage does not emit provenance (known gap)"


@dataclass
class Finding:
    level: Level
    code: str
    message: str
    fix: str = ""

    def as_dict(self) -> dict[str, str]:
        return {"level": self.level, "code": self.code, "message": self.message, "fix": self.fix}


@dataclass
class Stage:
    """One fixed-path stage of the pipeline.

    `inputs` and `outputs` are repo-relative. `provenance` lists outputs that
    must have a sibling `<stem>.provenance.json`; `provenance_exempt` lists
    outputs whose producing stage is known not to write one.
    """

    n: int
    name: str
    command: str
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    provenance: list[str] = field(default_factory=list)
    provenance_exempt: list[str] = field(default_factory=list)


# The fixed-path portion of the contract. Stages 1-3 (parse hypothesis, build
# test spec, derive search spec) write into a caller-chosen --outdir and so
# cannot be checked from a fixed path; they are reported separately.
CONTRACT: list[Stage] = [
    Stage(
        n=4,
        name="fetch ImmPort studies",
        command="hypothesis2omics_fetch_immport_studies(study_accessions=[...])",
        outputs=["data/immport_cache/manifest.json"],
    ),
    Stage(
        n=5,
        name="parse ImmPort Tab archives",
        command="hypothesis2omics_parse_immport_studies()",
        inputs=["data/immport_cache/manifest.json"],
        outputs=[
            "data/immport_cache/parsed/sample_manifest.tsv",
            "data/immport_cache/parsed/batch_manifest.json",
        ],
        provenance=["data/immport_cache/parsed/sample_manifest.tsv"],
    ),
    Stage(
        n=6,
        name="plan GEO retrieval",
        command="hypothesis2omics_plan_geo_retrieval()",
        inputs=["data/immport_cache/parsed/sample_manifest.tsv"],
        outputs=["data/geo_cache/plan/geo_download_plan.tsv"],
        provenance=["data/geo_cache/plan/geo_download_plan.tsv"],
    ),
    Stage(
        n=7,
        name="fetch planned GEO series",
        command="hypothesis2omics_fetch_planned_geo()",
        inputs=["data/geo_cache/plan/geo_download_plan.tsv"],
        outputs=["data/geo_cache/manifest.json"],
    ),
    Stage(
        n=8,
        name="parse GEO series matrices",
        command="hypothesis2omics_parse_geo_matrices()",
        inputs=["data/geo_cache/plan/geo_download_plan.tsv", "data/geo_cache/manifest.json"],
        outputs=["data/geo_cache/parsed/geo_matrix_parse_manifest.json"],
    ),
    Stage(
        n=9,
        name="build the validator-input bundle",
        command="python data/validator_handoff_parse_module.py data/validator_handoff_config.json",
        inputs=["data/geo_cache/parsed/geo_matrix_parse_manifest.json"],
        outputs=[
            "data/validator_input/sample_manifest.tsv",
            "data/validator_input/feature_expression.tsv",
            "data/validator_input/quantitative_outcome.tsv",
        ],
        provenance=[
            "data/validator_input/sample_manifest.tsv",
            "data/validator_input/feature_expression.tsv",
            "data/validator_input/quantitative_outcome.tsv",
        ],
    ),
    Stage(
        n=10,
        name="extract criterion-level evidence",
        command=(
            "python scientific_validator/handoff_adapter.py --input-dir data/validator_input "
            "--study-accession SDY1264 --gene EIF2AK4 "
            "--output scientific_validator/output/SDY1264_validator_input_evidence.json"
        ),
        inputs=["data/validator_input/sample_manifest.tsv"],
        outputs=["scientific_validator/output/SDY1264_validator_input_evidence.json"],
    ),
    Stage(
        n=11,
        name="assess scientific eligibility",
        command=(
            "python scientific_validator/eligibility_engine.py "
            "--spec scientific_validator/examples/yf17d_test_spec.json "
            "--dataset scientific_validator/output/SDY1264_validator_input_evidence.json"
        ),
        inputs=["scientific_validator/output/SDY1264_validator_input_evidence.json"],
        outputs=[],
    ),
]

# Bundles that are read directly by an analysis script rather than produced by a
# numbered stage. run_yf17d.py reads these, NOT data/validator_input/.
SIDE_BUNDLES: list[tuple[str, list[str]]] = [
    ("data/validator_immport", ["sample_manifest.tsv", "quantitative_outcome.tsv"]),
    ("data/validator_geo", ["feature_expression.tsv"]),
]

# Outputs whose producing script writes no provenance file. Recorded so the
# check does not cry wolf, and so the gap stays visible.
KNOWN_NO_PROVENANCE: list[str] = [
    "data/galaxy_file_input/*_intensities.csv",
]

# Directories every stage writes into. Scanned wholesale rather than by declared
# filename, because hand-rolled output is exactly the case that does NOT use a
# name the contract knows about.
SCAN_DIRS: list[str] = [
    "data/validator_input",
    "data/validator_immport",
    "data/validator_geo",
    "data/immport_cache/parsed",
    "data/geo_cache/plan",
]
SCAN_SUFFIXES = (".tsv", ".csv")


def run_check(cmd: list[str]) -> tuple[int, str]:
    """Run a short command, returning (returncode, combined output)."""
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=False)
    except FileNotFoundError:
        return 127, f"not found: {cmd[0]}"
    except subprocess.TimeoutExpired:
        return 124, "timed out"
    except OSError as err:
        return 1, str(err)
    return proc.returncode, (proc.stdout + proc.stderr).strip()


def check_environment(repo: Path) -> list[Finding]:
    """Can the stages actually run, and which tools are unavailable without credentials."""
    out: list[Finding] = []

    if shutil.which("uv"):
        out.append(Finding("ok", "uv", "uv is on PATH"))
    else:
        out.append(
            Finding(
                "block",
                "uv",
                "uv is not on PATH; every MCP server launches as `uv run`",
                "curl -LsSf https://astral.sh/uv/install.sh | sh",
            )
        )

    venv_python = repo / ".venv" / "bin" / "python"
    if not venv_python.exists():
        venv_python = repo / ".venv" / "Scripts" / "python.exe"
    if not venv_python.exists():
        out.append(
            Finding(
                "note",
                "venv",
                ".venv is absent; the scientific stack is ~230 MB and installing it on first "
                "launch outruns the agent's MCP startup window",
                f"uv sync --directory {repo}",
            )
        )
    else:
        # Touch a real attribute, not just the import. Both of these packages are
        # currently installed in a state where `import X` succeeds and the module is
        # a hollow namespace -- an import-only check reports them healthy.
        for module, importer, needed_by in (
            ("pandas", "import pandas; pandas.read_csv", "handoff_adapter.py, run_yf17d.py"),
            ("scipy.stats", "from scipy import stats; stats.pearsonr", "run_yf17d.py"),
        ):
            code, detail = run_check([str(venv_python), "-c", importer])
            if code == 0:
                out.append(Finding("ok", f"import:{module}", f".venv imports {module}"))
            else:
                out.append(
                    Finding(
                        "block",
                        f"import:{module}",
                        f".venv cannot use {module}: "
                        f"{detail.splitlines()[-1] if detail else code}; needed by {needed_by}",
                        f"uv sync --reinstall-package {module.split('.')[0]} --directory {repo}",
                    )
                )

    if os.environ.get("IMMPORT_API_KEY") or os.environ.get("IMMPORT_API_KEY_FILE"):
        out.append(Finding("ok", "immport-cred", "ImmPort credential present"))
    else:
        out.append(
            Finding(
                "note",
                "immport-cred",
                "no IMMPORT_API_KEY or IMMPORT_API_KEY_FILE; "
                "hypothesis2omics_fetch_immport_studies is unavailable, the rest of the "
                "pipeline is not affected",
                "export IMMPORT_API_KEY_FILE=data/immport_cache/<your-key>.json",
            )
        )

    if os.environ.get("OPENAI_API_KEY") or os.environ.get("ARGO_USER"):
        out.append(Finding("ok", "llm-cred", "LLM gateway credential present"))
    else:
        out.append(
            Finding(
                "note",
                "llm-cred",
                "no OPENAI_API_KEY or ARGO_USER; hypothesis_parser_parse_hypothesis and "
                "_build_test_spec are unavailable, every other stage is not affected",
                "export OPENAI_API_KEY=ac.yourname   # Argo username, or a real key",
            )
        )

    if os.environ.get("OPENAI_BASE_URL") and not os.environ.get("OPENAI_MODEL"):
        out.append(
            Finding(
                "block",
                "llm-model",
                "OPENAI_BASE_URL is set but OPENAI_MODEL is not, so the default model id "
                "'claudeopus5' applies -- that is an Argo id and means nothing elsewhere",
                "export OPENAI_MODEL=<a model the gateway at OPENAI_BASE_URL serves>",
            )
        )

    agent_dir = Path(os.environ.get("PI_CODING_AGENT_DIR") or (Path.home() / ".pi" / "agent"))
    mcp_json = agent_dir / "mcp.json"
    wanted = {"hypothesis-parser", "keyword-to-immport", "hypothesis2omics"}
    if not mcp_json.exists():
        out.append(
            Finding(
                "note",
                "mcp-registration",
                f"{mcp_json} does not exist, so Loom has no servers from this repo "
                "(irrelevant outside Loom)",
                "node mcp/register.mjs",
            )
        )
    else:
        try:
            servers = set(json.loads(mcp_json.read_text(encoding="utf-8")).get("mcpServers", {}))
        except (OSError, json.JSONDecodeError) as err:
            out.append(Finding("note", "mcp-registration", f"could not read {mcp_json}: {err}"))
        else:
            missing = sorted(wanted - servers)
            if missing:
                out.append(
                    Finding(
                        "note",
                        "mcp-registration",
                        f"{mcp_json} is missing {', '.join(missing)}",
                        "node mcp/register.mjs",
                    )
                )
            else:
                out.append(
                    Finding("ok", "mcp-registration", f"all three servers registered in {mcp_json}")
                )

    return out


def newest_mtime(repo: Path, rels: Iterable[str]) -> float:
    times = [(repo / r).stat().st_mtime for r in rels if (repo / r).exists()]
    return max(times) if times else 0.0


def check_state(repo: Path) -> tuple[list[Finding], list[Stage]]:
    """Report which artifacts exist, flag stale ones, and find every runnable stage."""
    out: list[Finding] = []
    runnable: list[Stage] = []

    for stage in CONTRACT:
        have_in = [r for r in stage.inputs if (repo / r).exists()]
        have_out = [r for r in stage.outputs if (repo / r).exists()]
        inputs_ready = len(have_in) == len(stage.inputs)
        complete = bool(stage.outputs) and len(have_out) == len(stage.outputs)

        if complete:
            in_t, out_t = newest_mtime(repo, stage.inputs), newest_mtime(repo, stage.outputs)
            if in_t and out_t and out_t < in_t:
                out.append(
                    Finding(
                        "note",
                        f"stale:{stage.n}",
                        f"stage {stage.n} ({stage.name}) output is older than its input",
                        stage.command,
                    )
                )
            else:
                out.append(Finding("ok", f"stage:{stage.n}", f"stage {stage.n} {stage.name}: done"))
        elif inputs_ready:
            detail = (
                "inputs ready, reports to stdout" if not stage.outputs else "inputs ready, not run"
            )
            out.append(
                Finding("note", f"stage:{stage.n}", f"stage {stage.n} ({stage.name}): {detail}",
                        stage.command)
            )
            runnable.append(stage)
        else:
            missing = sorted(set(stage.inputs) - set(have_in))
            out.append(
                Finding(
                    "note",
                    f"stage:{stage.n}",
                    f"stage {stage.n} ({stage.name}): blocked, missing {', '.join(missing)}",
                )
            )

    for directory, names in SIDE_BUNDLES:
        present = [n for n in names if (repo / directory / n).exists()]
        if present:
            out.append(
                Finding(
                    "ok",
                    f"bundle:{directory}",
                    f"{directory}/ holds {', '.join(present)} "
                    "(read by evidence_rules/run_yf17d.py, NOT data/validator_input/)",
                )
            )

    return out, runnable


def check_provenance(repo: Path) -> list[Finding]:
    """Flag stage outputs with no sibling provenance -- the signature of hand-rolled work."""
    out: list[Finding] = []
    checked = 0
    seen: set[str] = set()

    def verify(rel: str, producer: str) -> None:
        nonlocal checked
        target = repo / rel
        if not target.exists() or rel in seen:
            return
        seen.add(rel)
        checked += 1
        sibling = target.with_suffix("").with_suffix(".provenance.json")
        if sibling.suffix != ".json":
            sibling = target.parent / f"{target.stem}.provenance.json"
        if sibling.exists():
            return
        out.append(
            Finding(
                "block",
                "orphan",
                f"{rel} exists with no {sibling.relative_to(repo)}; "
                "a pipeline stage always writes one, so this artifact was produced outside "
                "the pipeline",
                f"regenerate it with: {producer}",
            )
        )

    for stage in CONTRACT:
        for rel in stage.provenance:
            verify(rel, stage.command)

    for directory, names in SIDE_BUNDLES:
        for name in names:
            verify(f"{directory}/{name}", "the stage that produced this bundle")

    for rel_dir in SCAN_DIRS:
        directory = repo / rel_dir
        if not directory.is_dir():
            continue
        for path in sorted(directory.iterdir()):
            if path.suffix in SCAN_SUFFIXES and path.is_file():
                verify(
                    str(path.relative_to(repo)).replace(os.sep, "/"),
                    "the stage that writes this directory -- see references/file-contract.md",
                )

    if checked and not out:
        out.append(
            Finding("ok", "orphan", f"{checked} stage outputs checked, all carry provenance")
        )
    for pattern in KNOWN_NO_PROVENANCE:
        if list(repo.glob(pattern)):
            out.append(
                Finding("note", "no-provenance-stage", f"{pattern}: {PROVENANCE_EXEMPT_NOTE}")
            )
    return out


SYMBOL = {"block": "FAIL", "note": "note", "ok": "ok  "}


def render(sections: list[tuple[str, list[Finding]]], runnable: list[Stage]) -> str:
    lines: list[str] = []
    for title, findings in sections:
        lines.append("")
        lines.append(title)
        lines.append("-" * len(title))
        for f in findings:
            lines.append(f"  [{SYMBOL[f.level]}] {f.message}")
            if f.fix and f.level != "ok":
                lines.append(f"           -> {f.fix}")
    lines.append("")
    if runnable:
        lines.append("Runnable now:")
        for stage in runnable:
            lines.append(f"  {stage.n} -- {stage.name}")
            lines.append(f"       {stage.command}")
    else:
        lines.append("Runnable now: nothing; every fixed-path stage has its outputs.")
    lines.append("")
    lines.append(
        "Stages 1-3 (parse hypothesis, build test spec, derive search spec) write into a "
        "caller-chosen outdir and are not checked here."
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--repo", type=Path, default=Path.cwd(), help="checkout of hypothesis2omics"
    )
    parser.add_argument("--json", action="store_true", help="emit JSON instead of text")
    args = parser.parse_args(argv)

    repo: Path = args.repo.resolve()
    if not (repo / "pyproject.toml").exists() or not (repo / "evidence_rules").is_dir():
        print(f"error: {repo} does not look like a hypothesis2omics checkout", file=sys.stderr)
        return 2

    env = check_environment(repo)
    state, runnable = check_state(repo)
    prov = check_provenance(repo)
    sections = [
        ("Environment", env),
        ("Pipeline state", state),
        ("Provenance", prov),
    ]
    blocking = [f for _, fs in sections for f in fs if f.level == "block"]

    if args.json:
        payload: dict[str, Any] = {
            "repo": str(repo),
            "exit": 1 if blocking else 0,
            "runnable": [
                {"n": s.n, "name": s.name, "command": s.command} for s in runnable
            ],
            "sections": {title: [f.as_dict() for f in fs] for title, fs in sections},
        }
        print(json.dumps(payload, indent=2))
    else:
        print(render(sections, runnable))

    return 1 if blocking else 0


if __name__ == "__main__":
    raise SystemExit(main())
