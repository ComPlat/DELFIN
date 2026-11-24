"""What the request asked of the model, and whether anyone asked it.

Input: the session's settings and environment. Output: the sampling
parameters to put in a request, and a record of what was actually sent.

THE DEFECT THIS CLOSES. The run manifest records ``model.temperature``
and ``model.seed``, and a survey of 26 autonomous-research systems
(arXiv:2608.05179) counts released seeds and traces as the field's
largest gap -- 38% against 83% for code. Our manifest asked the engine
for both, the engine has neither, and ``temperature`` and ``seed``
appear NOWHERE in this product: nothing here has ever set them. So both
fields were ``unresolved`` in every manifest, which reads as "we forgot
to record it" when the truth is "nobody asked for it and the provider's
default applied".

Those are different facts and the second one is worth recording. A run
at an unstated default temperature is NOT reproducible, and saying so is
more useful than an empty field that looks like an oversight.

WHAT THIS DOES NOT DO. It does not make runs deterministic by default.
Sending ``temperature: 0`` to every endpoint would change what every
existing run does, and that decision belongs to whoever is running the
experiment, not to a recording module. Default: nothing is sent, the
provider's default applies, and the record says exactly that.

A benchmark that wants reproducible sampling sets it -- which is what a
pre-registration should specify -- and the manifest then carries the
values with ``source: "configured"``.

AND IT MUST NOT COST FUNCTION. An endpoint that rejects a parameter it
does not know answers 400, and a long unattended run must not die of a
recording feature. ``without(params)`` gives the caller the retry
payload, and ``REFUSED`` is the source recorded afterwards, so the run
continues and the manifest says the endpoint refused the parameter
rather than claiming it was applied.
"""
from __future__ import annotations

import os
from typing import Any, Optional

#: Nothing was asked: the provider's own default applied. Not an
#: oversight, and not reproducible either.
PROVIDER_DEFAULT = "provider_default"
#: The session asked for these values and the request carried them.
CONFIGURED = "configured"
#: The session asked and the endpoint rejected the parameter (400). The
#: values are reported as what was asked, never as what applied.
REFUSED = "refused_by_endpoint"

#: Environment overrides, so one benchmark run can fix sampling without
#: touching the user's settings file.
ENV_TEMPERATURE = "DELFIN_SAMPLING_TEMPERATURE"
ENV_SEED = "DELFIN_SAMPLING_SEED"


#: Set once an endpoint has refused the parameters, so the rest of the
#: process stops sending them. Per-process and deliberately not
#: persisted: the next session may run against another endpoint, and a
#: refusal remembered across sessions would silently stop asking.
_disabled = False


def disable() -> None:
    """Stop offering the parameters for the rest of this process."""
    global _disabled
    _disabled = True


def disabled() -> bool:
    return _disabled


def _reset_for_tests() -> None:
    """Clear the process flag. Tests only -- a test that disabled
    sampling must not decide the next test's answer."""
    global _disabled
    _disabled = False


def _as_float(raw: Any) -> Optional[float]:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    # A temperature outside the range every endpoint accepts is a
    # configuration error, and sending it would be a 400 on every turn.
    return value if 0.0 <= value <= 2.0 else None


def _as_int(raw: Any) -> Optional[int]:
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def configured(settings: Optional[dict] = None) -> dict[str, Any]:
    """The sampling the session asks for, or empty for "nothing".

    The environment wins over the settings file: a benchmark arm fixes
    sampling for its own process without editing the user's settings.
    """
    if _disabled:
        # An endpoint already refused them; asking again costs a round
        # per turn and answers the same way.
        return {}
    out: dict[str, Any] = {}
    agent = {}
    if settings is None:
        try:
            from delfin import user_settings
            agent = (user_settings.load_settings() or {}).get("agent") or {}
        except Exception:                               # noqa: BLE001
            agent = {}
    else:
        agent = settings.get("agent") or settings

    temperature = _as_float(os.environ.get(ENV_TEMPERATURE))
    if temperature is None:
        temperature = _as_float(agent.get("sampling_temperature"))
    if temperature is not None:
        out["temperature"] = temperature

    seed = _as_int(os.environ.get(ENV_SEED))
    if seed is None:
        seed = _as_int(agent.get("sampling_seed"))
    if seed is not None:
        out["seed"] = seed
    return out


def record(params: dict[str, Any], *, source: str = "") -> dict[str, Any]:
    """What the manifest should carry about this request's sampling.

    Every field is present and named, so a reader never has to decide
    whether an absent key means "default" or "forgotten".
    """
    if not source:
        source = CONFIGURED if params else PROVIDER_DEFAULT
    return {
        "temperature": params.get("temperature"),
        "seed": params.get("seed"),
        "source": source,
        "reproducible": bool(params.get("seed") is not None
                             and params.get("temperature") == 0.0
                             and source == CONFIGURED),
    }


def without(payload: dict[str, Any]) -> dict[str, Any]:
    """*payload* with the sampling parameters removed, for the retry.

    An endpoint that does not know ``seed`` answers 400 to the whole
    request. A run must not end because of a field added for the record,
    so the caller retries without them and records REFUSED.
    """
    out = dict(payload)
    out.pop("temperature", None)
    out.pop("seed", None)
    return out


def looks_like_a_parameter_refusal(exc: BaseException) -> bool:
    """Whether *exc* is an endpoint rejecting one of these parameters.

    Deliberately narrow: it must mention one of the parameters by name.
    A 400 about the messages, the model or the context is a real error
    and retrying without sampling would hide it.
    """
    text = f"{type(exc).__name__}: {exc}".lower()
    if "400" not in text and "invalid" not in text and "unsupported" not in text:
        return False
    return "temperature" in text or "seed" in text
