"""The provider guessed from a base URL follows its host, not a substring.

A substring test sent any URL that merely contained ``kit.edu`` (in a path
or query) the KIT credentials. The host decides; the legitimate URLs keep
their old answers.
"""

import pytest

from delfin.agent.api_client import _infer_provider_from_base_url as infer


@pytest.mark.parametrize("url, want", [
    ("https://ki-toolbox.scc.kit.edu/api/v1", "kit"),
    ("https://KI-TOOLBOX.SCC.KIT.EDU/api/v1", "kit"),
    ("http://localhost:11434/v1", "ollama"),
    ("http://127.0.0.1:11434", "ollama"),
    ("http://localhost:8000/v1", "ollama"),
    ("localhost:11434", "ollama"),
    ("https://api.openai.com/v1", "openai"),
    ("", "openai"),
])
def test_known_urls_keep_their_provider(url, want):
    assert infer(url) == want


@pytest.mark.parametrize("url", [
    "https://other.example/?kit.edu",
    "https://other.example/kit.edu/v1",
    "https://kit.edu.other.example/v1",
    "https://other.example/ki-toolbox",
    "https://notkit.edu.example.org",
])
def test_kit_in_path_query_or_foreign_host_is_not_kit(url):
    assert infer(url) == "openai"
