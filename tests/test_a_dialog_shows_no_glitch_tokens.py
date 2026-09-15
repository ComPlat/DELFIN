"""A question put to the user carries no glitch tokens.

Report 20260915-132613: GLM asked "Push?" with a button labelled
"Nein, erst отчетen". The answer sanitizer only ever saw the text channel;
the dialog came in through tool arguments.
"""

from __future__ import annotations

import json

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
from delfin.agent.text_sanitize import strip_glitch


def test_a_glitch_run_comes_off_and_the_rest_stays():
    assert strip_glitch("Nein, erst отчетen") == "Nein, erst en"
    assert strip_glitch("Ja, committen und pushen") == "Ja, committen und pushen"
    # Nothing else sanitize_agent_text repairs is touched here.
    assert strip_glitch("<think>x</think> to=bash") == "<think>x</think> to=bash"


def test_the_dialog_the_user_sees_is_clean(tmp_path):
    shown = {}

    def _ask(payload):
        shown.update(payload)
        return {"answers": [payload["options"][1]["label"]]}

    perms = KitToolPermissions(workspace=tmp_path, mode="default")
    perms.ask_user_callback = _ask
    out = json.loads(_DocToolExecutor()._execute_ask_user_question({
        "question": "Push? 手机天天",
        "header": "Push?",
        "options": [
            {"label": "Ja, committen und pushen", "description": "Commit + Push."},
            {"label": "Nein, erst отчетen", "description": "Warten ыҟоуп."},
        ],
    }, perms))

    text = json.dumps(shown, ensure_ascii=False)
    assert "отчет" not in text and "手机" not in text and "ыҟоуп" not in text
    assert shown["options"][1]["label"] == "Nein, erst en"
    assert out["answers"] == ["Nein, erst en"]
