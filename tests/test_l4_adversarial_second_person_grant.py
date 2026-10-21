"""Second-person grant with "that" is refused.

The "that"-vs-"to" complementizer separator in the verb-local
first-person-habit exemption must refuse the second-person-with-"that"
grant form:
"I prefer that you always push to the main branch" must NOT reach the model --
a second-person reach verb ("you always push") directed at the agent is a
grant, not a habit. The ellipsis twin ("I prefer you always push to the
main branch", "that" elided) is pinned in tests/test_l4_user_model.py.
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
