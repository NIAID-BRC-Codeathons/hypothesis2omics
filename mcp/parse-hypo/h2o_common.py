#!/usr/bin/env python3
"""Shared machinery for the hypothesis -> spec pipeline.

Everything here is repository-agnostic: it talks to an OpenAI-compatible gateway
and renders YAML, and it imports nothing from the ImmPort layer. That is deliberate
-- it lets parse_hypothesis.py run and be tested with no ImmPort dependency at all.

    h2o_common.py          <- this file: LLM calls, JSON extraction, YAML scalars
        ^         ^
        |         |
    parse_hypothesis.py    step 1, offline-testable
        ^
        |
    build_test_spec.py     steps 2-3, imports the ImmPort layer

Install:
    pip install openai httpx pyyaml

Not a script; imported by the two step modules.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from typing import Any

# Any OpenAI-compatible gateway. The default is Argo, so an existing ARGO_USER
# setup keeps working untouched; point OPENAI_BASE_URL at something else and this
# is an ordinary OpenAI client. OPENAI_URL is accepted as an alias for it.
ARGO_BASE_URL = "https://apps.inside.anl.gov/argoapi/v1"
OPENAI_BASE_URL = (
    os.environ.get("OPENAI_BASE_URL") or os.environ.get("OPENAI_URL") or ARGO_BASE_URL
)
# On Argo this is a username rather than a secret, which is why ARGO_USER still
# works. Either way it lands in shell history, so keep it in a per-project file
# you source, not in ~/.bashrc or a committed config.
OPENAI_API_KEY = (
    os.environ.get("OPENAI_API_KEY") or os.environ.get("ARGO_USER") or "ac.jdoe"
)
# `or`, never a get() default: an MCP client that forwards an unset variable hands
# us "" rather than leaving it absent, and "" is a value get() would hand straight
# back. Every lookup here treats empty as missing for that reason.
DEFAULT_MODEL = os.environ.get("OPENAI_MODEL") or "claudeopus5"
# claudeopus5 returns empty completions on some ordinary biomedical inputs (see
# _complete). Anything it drops is retried here, and provenance records which
# model actually produced the artifact. That quirk is Argo's, and "claudesonnet5"
# is not a model id anywhere else, so the retry turns itself off as soon as you
# name your own model. Set OPENAI_FALLBACK_MODEL to pick a different second try.
FALLBACK_MODEL = os.environ.get("OPENAI_FALLBACK_MODEL") or (
    "claudesonnet5" if not os.environ.get("OPENAI_MODEL") else None
)
# Ceiling for the automatic budget growth in _complete_json. Kept under the
# 21000 cap Argo imposes on streaming Anthropic requests, so the same number
# stays valid if this ever streams.
MAX_TOKEN_BUDGET = 16000
GENERATOR_VERSION = "0.4.0"

YAML_BOOLEANS = frozenset("y n yes no true false on off".split())


def client() -> Any:
    """An OpenAI client for the configured gateway. Lazy, so --help works offline."""
    from openai import OpenAI

    return OpenAI(api_key=OPENAI_API_KEY, base_url=OPENAI_BASE_URL)


# ------------------------------------------------------------------------ LLM


def extract_json(text: str) -> dict[str, Any]:
    """Pull one JSON object out of a model response, fenced or not."""
    text = text.strip()
    fence = re.search(r"```(?:json)?\s*(.*?)```", text, re.DOTALL)
    if fence:
        text = fence.group(1).strip()
    start, end = text.find("{"), text.rfind("}")
    if start == -1 or end <= start:
        raise ValueError(f"no JSON object in model response: {text[:200]!r}")
    return json.loads(text[start : end + 1])


class EmptyCompletion(RuntimeError):
    """The gateway returned a successful response with no content.

    `completion_tokens` separates two failure modes that need opposite fixes:
    zero means the model generated nothing at all and a different model is the
    only way forward; non-zero means it generated something internally and the
    budget ran out, so a larger max_tokens is the fix.
    """

    def __init__(self, message: str, completion_tokens: int | None = None):
        super().__init__(message)
        self.completion_tokens = completion_tokens


def _complete(client, model: str, messages: list[dict], max_tokens: int) -> str:
    resp = client.chat.completions.create(
        model=model,
        messages=messages,
        max_tokens=max_tokens,
        # Argo rejects every value other than 1 for claudesonnet5 and the gpt5
        # family, and the SDK injects a default when the argument is absent.
        temperature=1,
    )
    choice = resp.choices[0]
    content = choice.message.content or ""
    if content.strip():
        return content

    tokens = getattr(resp.usage, "completion_tokens", None)
    if tokens:
        # Text budget exhausted by internal reasoning. Measured 2026-09-16:
        # claudesonnet5 returned finish_reason 'stop' with completion_tokens
        # equal to the full max_tokens and no content. The same request succeeds
        # with a larger budget, so this is not a model-level failure.
        raise EmptyCompletion(
            f"{model} produced no text but consumed {tokens} completion tokens "
            f"(finish_reason={choice.finish_reason!r}). The budget went to internal "
            f"reasoning; retry with a larger max_tokens.",
            tokens,
        )

    # Measured 2026-09-16: claudeopus5 on Argo returns an empty completion --
    # finish_reason 'stop', completion_tokens 0, no error and no refusal text --
    # for some ordinary biomedical inputs. "Depletion of CD4+ T cells causes a
    # reduction in germinal center B-cell formation after immunization." was
    # empty 8/8 while a control prompt was empty 0/8; claudesonnet5 answered the
    # same input normally. The cause is not observable from outside the gateway.
    raise EmptyCompletion(
        f"{model} returned an empty completion with zero completion tokens "
        f"(finish_reason={choice.finish_reason!r}). This is input-specific and "
        f"reproducible, not transient, and a larger budget will not help. "
        f"A different model is the only fix.",
        0,
    )


def _complete_json(client, model: str, messages: list[dict], max_tokens: int,
                   attempts: int = 2) -> dict:
    """_complete plus JSON extraction, retried once on unparseable output.

    Measured 2026-09-16: with temperature forced to 1, responses occasionally
    carry an unescaped double quote inside a string value, which is invalid JSON.
    It is sampling noise rather than a prompt defect -- the identical request
    parses on a retry -- so one retry is cheaper than a tolerant parser that
    might silently mis-read a genuinely malformed object.
    """
    last: Exception | None = None
    budget = max_tokens
    for _ in range(attempts):
        try:
            return extract_json(_complete(client, model, messages, budget))
        except EmptyCompletion as exc:
            # Zero tokens is a model-level dead end; let the caller fall back to
            # another model. Any other count means the budget was the problem.
            if not exc.completion_tokens:
                raise
            budget = min(budget * 2, MAX_TOKEN_BUDGET)
            print(f"note: {exc} -- retrying with max_tokens={budget}", file=sys.stderr)
            last = exc
        except ValueError as exc:  # JSONDecodeError subclasses ValueError
            # A truncated response is indistinguishable from a malformed one at
            # this level, and both are cheap to retry, so grow the budget too.
            budget = min(budget * 2, MAX_TOKEN_BUDGET)
            last = exc
    raise ValueError(f"{model} gave no usable JSON in {attempts} attempts: {last}")


def _ask(client, model, system, user, validator, max_tokens, label,
         fallback: str | None = FALLBACK_MODEL) -> tuple[dict, dict]:
    """One call plus, if the validator objects, one repair round-trip.

    Returns the object and a dict of what it took to get it: which model
    actually answered, and whether a repair round-trip was needed.
    """
    messages = [{"role": "system", "content": system}, {"role": "user", "content": user}]
    used = model
    try:
        obj = _complete_json(client, model, messages, max_tokens)
    except EmptyCompletion as exc:
        if not fallback or fallback == model:
            raise
        print(f"note: {exc}", file=sys.stderr)
        print(f"note: {label} falling back to {fallback}", file=sys.stderr)
        used = fallback
        obj = _complete_json(client, used, messages, max_tokens)

    problems = validator(obj)
    meta = {"model_used": used, "repair_pass_used": False}
    if not problems:
        return obj, meta

    print(f"note: {label} had {len(problems)} problem(s); requesting repair", file=sys.stderr)
    messages += [
        # The re-serialised object, not the raw text: it is what the validator
        # actually judged, and it is guaranteed to be valid JSON.
        {"role": "assistant", "content": json.dumps(obj)},
        {
            "role": "user",
            "content": "That JSON was rejected by the downstream validator:\n"
            + "\n".join(f"- {p}" for p in problems)
            + "\n\nReturn the corrected JSON object, nothing else.",
        },
    ]
    obj = _complete_json(client, used, messages, max_tokens)
    problems = validator(obj)
    if problems:
        raise ValueError(f"{label} still invalid after repair:\n" + "\n".join(problems))
    meta["repair_pass_used"] = True
    return obj, meta


# ----------------------------------------------------------------------- YAML


def _scalar(value: Any) -> str:
    """Render a scalar as YAML, quoting only when a bare word would misparse."""
    if value is None:
        return "null"
    text = str(value)
    safe = re.fullmatch(r"[A-Za-z][A-Za-z0-9 _.+/()'\[\]-]*", text)
    if safe and text.lower() not in YAML_BOOLEANS and "  " not in text:
        return text
    return json.dumps(text)


def _header(provenance: dict[str, Any] | None) -> list[str]:
    return [f"# {k}: {v}" for k, v in provenance.items()] + [""] if provenance else []


# ----------------------------------------------------------------- provenance


def _provenance(generator: str, requested: str, system: str, step: str,
                meta: dict) -> dict[str, Any]:
    prov = {
        "step": step,
        "generator": f"{generator} {GENERATOR_VERSION}",
        "model": meta["model_used"],
        "temperature": 1,
        "prompt_sha256": hashlib.sha256(system.encode()).hexdigest()[:16],
        "repair_pass_used": meta["repair_pass_used"],
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    if meta["model_used"] != requested:
        # The artifact must never claim a model that did not produce it.
        prov["model_requested"] = requested
        prov["fallback_reason"] = "empty completion"
    return prov
