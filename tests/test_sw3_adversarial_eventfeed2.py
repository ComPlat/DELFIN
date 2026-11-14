"""Adversarial tests for LeadEventFeed throttling and payload robustness.

Reviewer nacht-s21, package SW3. Attacks on delfin/agent/swarm/lead.py:

1. An event dropped because the throttle window is full is marked seen
   BEFORE the throttle decision (lead.py:253 before lead.py:258) -- so a
   retry of that event after the window resets is swallowed forever.
   The dedup is documented (lead.py:196-204) as catching retries of the
   SAME event, but a throttle drop is not a retry.
2. _event_identity serialises the payload with json.dumps(...,
   sort_keys=True) (lead.py:211). HYPOTHESIS REFUTED by measurement:
   json.dumps coerces non-string keys without raising, so a payload
   with an int key flows through and is delivered like any other.
   The test pins that behaviour instead (no crash, later events flow).
"""
from delfin.agent.swarm import lead as lead_mod


class _Clock:
    def __init__(self) -> None:
        self.value = 1000.0

    def __call__(self) -> float:
        return self.value


def _event(n: int) -> dict:
    return {"kind": "phase_change", "package": "PK", "member": f"m{n}",
            "phase": n}


def test_event_dropped_by_throttle_is_delivered_when_reemitted():
    """A distinct event suppressed by the throttle must not be swallowed
    forever: when the emitter re-sends it after the window resets, it
    must be delivered (it is not a retry of a delivered event)."""
    clock = _Clock()
    delivered = []
    feed = lead_mod.LeadEventFeed(max_per_window=2, window_s=10.0,
                                  now=clock,
                                  on_deliver=delivered.append)
    feed._on_event(_event(1))
    feed._on_event(_event(2))
    clock.value += 11.0  # window resets on the next event
    feed._on_event(_event(3))  # delivered: window empty again
    assert len(delivered) == 3, delivered
    # Now overflow the window: event 4 and 5 exceed max_per_window=2.
    feed._on_event(_event(4))
    feed._on_event(_event(5))
    assert len(delivered) == 4, delivered
    clock.value += 11.0
    # The emitter retries event 5 after the window reset. A throttle
    # drop is not a retry of a delivered event: the dedup key must NOT
    # have been recorded, so the retry must be delivered.
    feed._on_event(_event(5))
    assert len(delivered) == 5, delivered


def test_non_string_payload_keys_do_not_kill_the_feed():
    """Pinned behaviour (hypothesis 'crashes' refuted by measurement):
    json.dumps coerces non-string keys, so a payload with an int key is
    delivered like any other event, raises nothing, and later healthy
    events still flow."""
    delivered = []
    feed = lead_mod.LeadEventFeed(now=_Clock(), on_deliver=delivered.append)
    bad = {"kind": "ci_red", "package": "PK", "member": "m1",
           1: "non-string key"}
    feed._on_event(bad)  # must not raise
    feed._on_event(_event(7))
    assert len(delivered) == 2, delivered
