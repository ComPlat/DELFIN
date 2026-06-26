"""The flapping-KIT-model error (ServiceUnavailable) must NOT be mistaken for
an auth/401 problem — it echoes 'Received Model Group=…' which used to trip the
auth branch and tell the user their key was invalid (bug 2026-06-26)."""
from delfin.dashboard.tab_agent import _classify_model_error_text as C


def test_kit_service_unavailable_is_temp_not_auth():
    dump = ("litellm.ServiceUnavailableError: Hosted_vllmException - "
            "{\"error\":\"Model 'qwen3.5-397b-A17b' is temporarily unavailable. "
            "Please try again later.\"}. Received Model Group=kit.qwen3.5-397b-A17b")
    assert C(dump) == "temp"          # NOT "auth"


def test_real_auth_errors_are_auth():
    assert C("litellm.AuthenticationError: invalid api key") == "auth"
    assert C("Access denied for model group X") == "auth"
    assert C("HTTP 401: model not provisioned") == "auth"


def test_unrelated_error_is_neither():
    assert C("connection reset by peer") == ""
    assert C("") == ""
