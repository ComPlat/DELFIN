"""A gateway that lost its model server is not a bad request.

Measured 2026-09-21 at the KIT endpoint: a turn died on

    Error code: 400 - {'detail': 'Open WebUI: Server Connection Error'}

and the session stood at its prompt for 24 minutes. The SDK wraps it
as a BadRequestError and 400 is not a transient status, so nothing
retried it. The gateway reporting that ITS OWN connection to the model
server failed cannot be a malformed chat request -- the same shape as
the infrastructure markers beside it.

"Stopped while waiting for the endpoint", proposed alongside, is NOT
here on purpose: it is the notice DELFIN prints when someone stopped
the turn, and a retry would run the request the user just ended.
"""

from __future__ import annotations

from delfin.agent.api_client import _is_transient_api_error


class _BadRequest(Exception):
    status_code = 400


def test_the_measured_400_is_transient():
    assert _is_transient_api_error(_BadRequest(
        "Error code: 400 - {'detail': 'Open WebUI: Server Connection "
        "Error'}"))


def test_a_genuine_bad_request_still_fails_at_once():
    for msg in ("Error code: 400 - {'error': {'message': "
                "'The model `no-such-model` does not exist'}}",
                "Error code: 400 - {'error': {'message': "
                "'This model supports at most 8192 tokens'}}"):
        assert not _is_transient_api_error(_BadRequest(msg)), msg


# Measured 2026-09-25 on kit.glm-5.3: four supervised sessions ended
# their turns on these, each needing the operator to restart them.
_MEASURED_2026_09_25 = (
    "Error code: 400 - {'detail': \"litellm.APIConnectionError: "
    "APIConnectionError: Hosted_vllmException - Response payload is not "
    "completed: <TransferEncodingError: 400, message='Not enough data to "
    "satisfy transfer length header.'>\"}",
    "Error code: 400 - {'detail': 'No deployments available for selected "
    "model, Try again in 5 seconds. Passed model=kit.glm-5.3'}",
)


def test_a_response_cut_off_by_the_gateway_is_transient():
    for msg in _MEASURED_2026_09_25:
        assert _is_transient_api_error(_BadRequest(msg)), msg


def test_a_model_the_gateway_does_not_know_is_not_retried():
    # Same day, same outage -- but "model not found" is also what a typo
    # in the model name says, so it keeps failing at once.
    assert not _is_transient_api_error(_BadRequest(
        "Error code: 400 - {'detail': 'litellm.NotFoundError: NotFoundError: "
        "Hosted_vllmException - {\"error\":\"Model \\'glm-5.3\\' not found."
        "\"}'}"))
