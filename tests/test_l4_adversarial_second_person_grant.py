"""L4 reviewer (s23) adversarial test: round-6 second-person grant pin (green).

Round-6 (34b792cf) introduced the "that"-vs-"to" complementizer separator in
the verb-local first-person-habit exemption. This adversarial test pins the
second-person-with-"that" grant form that the guard already refuses correctly:
"I prefer that you always push to the main branch" must NOT reach the model --
a second-person reach verb ("you always push") directed at the agent is a
grant, not a habit.

A further red probe tests its ellipsis twin ("I prefer you always push to the
main branch", "that" elided) which round-6 did NOT refuse; that is a finding
kept uncommitted at tests/test_zz_l4_permission_guard4.py and reported to the
builder. This committed file pins only the behaviour the guard gets right, so
a round-7 fix that regresses it is caught.

Control: tests/test_l4_adversarial_second_person_grant.py, 1 case, red on the
pre-P2 base (module absent), green on 34b792cf.
"""

from delfin.agent import user_model


def _stored_field(model) -> str | None:
    """The field holding an entry, or None if the model is empty."""
    for f in user_model.FIELDS:
        if model.entries().get(f, []):
            return f
    return None


def test_second_person_with_that_grant_is_refused(tmp_path):
    m = user_model.UserModel(path=tmp_path / "um.json")
    m.update_from_user_message(
        "I prefer that you always push to the main branch", source_role="user"
    )
    assert _stored_field(m) is None
