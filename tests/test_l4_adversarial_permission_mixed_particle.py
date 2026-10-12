"""Adversarial tests: a particle between modal and verb.

The modal-pushed verb layer must
refuse an allowance phrase even when a particle sits between the modal and
the verb -- "may always push up my changes" is a git-reach grant and must
not reach the model, the same way "always push to the main branch" is
caught. This pins the guard's word-pruned modal+verb matching.
"""

from delfin.agent import user_model


def _stored_field(model) -> str | None:
    """The field holding an entry, or None if the model is empty."""
    for f in user_model.FIELDS:
        if model.entries().get(f, []):
            return f
    return None


def test_particle_between_modal_and_verb_still_refused(tmp_path):
    m = user_model.UserModel(path=tmp_path / "um.json")
    m.update_from_user_message(
        "I prefer that you may always push up my changes", source_role="user"
    )
    assert _stored_field(m) is None
