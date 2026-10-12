"""Adversarial test: a "you" after the reach clause.

The subject test scans the FULL left context before the frequency modal
for a _SECOND_PERSON token, which closes the compound-subject grants. This pins that the scan did NOT become
over-broad: a "you" that appears AFTER the reach clause (a separate addressed
clause, not the modal's subject) must not turn a genuine first-person habit
into a grant. "I prefer to always write my tests before you review them" is a
working-style habit and must store.

A post-modal "you" is exactly the class a future over-tightening of the
left-context scan could regress, so this is a canary.
"""

from delfin.agent import user_model


def _stored_field(model) -> str | None:
    """The field holding an entry, or None if the model is empty."""
    for f in user_model.FIELDS:
        if model.entries().get(f, []):
            return f
    return None


def test_post_modal_you_does_not_turn_habit_into_grant(tmp_path):
    m = user_model.UserModel(path=tmp_path / "um.json")
    m.update_from_user_message(
        "I prefer to always write my tests before you review them",
        source_role="user",
    )
    assert _stored_field(m) is not None
