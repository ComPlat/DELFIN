"""The check a lesson passes before it leaves the personal layer.

Input of :func:`check`: the lesson text, whether the target repository
is public, the identity of the running account (home, account, host;
default :func:`api_client._identity_for_commit_texts`), the names of the
project the lesson came from, the spans the person marked releasable in
the share dialog, and terms the person flagged as unpublished. Output: a
list of :class:`Finding` -- rule, reason, the offending snippet and its
offsets in the ORIGINAL text. Empty means the text may be offered.

Semantics: findings are reported, never repaired. A rewritten text is a
text the person did not read; the person rewords and the check runs
again on the new text. Any doubt blocks ("better ask more often than let
sensitive data go public").

Two rule sets. BASE applies to every repository layer: credentials
(``output_guard.scrub_secrets``), emails, the account's own home path,
account and host names (``output_guard.home_paths_in``, the core of the
commit-message identity check), machine references (absolute paths,
user@host, job ids: ``skill_share._machine_ref_snippets``), control
characters, licensed or internal markers (CSD/CCDC, refcode-shaped
tokens, "unpublished", "confidential", ...), user-specific phrasing
(first person, "yesterday", "this session"), internal URLs, and every
instruction about permissions, checks or protected files
(``skill_safety._TEXT_RULES`` on the raw and the folded text,
``untrusted.flags``, plus the check-disable and protected-zone rules
below). PUBLIC adds: generic host names (FQDNs, cluster node names),
@handles, decimal values and values with units, sample ids, every URL
outside a short allow-list, and the source project's names.

Releasable: only findings of the kinds in :data:`RELEASABLE_RULES`
(values, sample ids, refcode-shaped tokens, allow-list URLs) can be
released, by naming the exact snippet. A secret, an email, a home path
or an instruction cannot.

Rejected alternatives: redacting and passing (the person would approve
text that is not what gets published); a model-based classifier (not
deterministic, not testable without a network, and steerable by the very
text it judges).
"""
from __future__ import annotations

import difflib
import re
from dataclasses import dataclass
from typing import Iterable, Optional

__all__ = ["Finding", "check", "instruction_findings", "render",
           "RELEASABLE_RULES"]


@dataclass(frozen=True)
class Finding:
    rule: str
    reason: str
    snippet: str
    start: int = -1
    end: int = -1

    @property
    def releasable(self) -> bool:
        return self.rule in RELEASABLE_RULES

    def line(self) -> str:
        where = f" at {self.start}-{self.end}" if self.start >= 0 else ""
        return f"{self.rule}: {self.reason} -- {self.snippet!r}{where}"


#: Kinds the person may release in the dialog by naming the snippet.
RELEASABLE_RULES = frozenset({"number", "sample-id", "refcode-like", "url"})

# ---------------------------------------------------------------------------
# patterns
# ---------------------------------------------------------------------------

_EMAIL_RE = re.compile(r"(?<![\w.+-])[\w.+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)+")
_URL_RE = re.compile(r"\b[a-z][a-z0-9+.-]*://[^\s)>\]'\"]+", re.I)
#: Hosts a public lesson may cite: persistent identifiers only.
_PUBLIC_URL_HOSTS = ("doi.org", "arxiv.org")
_INTERNAL_HOST_RE = re.compile(
    r"^(?:localhost|127\.\d+\.\d+\.\d+|10\.\d+\.\d+\.\d+|192\.168\.\d+\.\d+"
    r"|172\.(?:1[6-9]|2\d|3[01])\.\d+\.\d+|\d+\.\d+\.\d+\.\d+"
    r"|[\w.-]+\.(?:local|internal|lan|corp|intra|intranet|home|test))$", re.I)

#: A licence or confidentiality marker. German stays German: it matches
#: what a German-writing user types.
_LICENSED_RE = re.compile(
    r"(?i)\b(?:CSD|CCDC|refcodes?|licen[sc]ed|proprietary|confidential"
    r"|internal[\s-]+only|unpublished|not\s+for\s+(?:distribution|release)"
    r"|under\s+NDA|NDA|embargo(?:ed)?|vertraulich|unver(?:ö|oe)ffentlicht"
    r"|nur\s+intern)\b")
#: Six capitals and an optional two-digit suffix: the CSD refcode shape.
_REFCODE_RE = re.compile(r"\b[A-Z]{6}(?:\d{2})?\b")
#: Upper-case words of that shape that are not refcodes.
_REFCODE_WORDS = frozenset({
    "DELFIN", "PYTHON", "GITHUB", "README", "CONFIG", "STDOUT", "STDERR",
    "SIGINT", "SIGHUP", "BASH", "UNICODE", "SCRIPT", "ERRORS", "OUTPUT",
    "STATUS", "SOURCE", "FAILED", "PASSED", "SLURMD", "NUMBER", "STRING",
    "OPTION", "BUFFER", "SERIAL", "SINGLE", "DOUBLE", "TRIPLE", "MANTA",
    "GOAT", "CREST", "DOCKER", "LINUX", "MACOS", "WINDOWS", "PYTEST",
    "SCRATCH", "PUBLIC", "PRIVATE", "IMPORT", "RETURN", "SELECT", "INSERT",
    "UPDATE", "DELETE", "CREATE", "BRANCH", "COMMIT", "MERGED", "REVIEW",
    "ALWAYS", "NEVER", "BEFORE", "AFTER", "SHOULD", "CANNOT",
})

#: First person, the reader addressed as the one who said something, and
#: time words that only mean something inside one session.
_USER_SPECIFIC_RE = re.compile(
    r"(?i)\b(?:I|I'm|I've|me|my|mine|we|we're|our|ours|you\s+(?:said|asked"
    r"|told|wanted)|the\s+user\s+(?:said|asked|wants|told)|yesterday|today"
    r"|tonight|this\s+morning|just\s+now|this\s+session|last\s+session"
    r"|ich|mein(?:e[mnrs]?)?|wir|unser(?:e[mnrs]?)?|heute|gestern"
    r"|vorhin)\b")

#: Switching off a check, in words. Complements skill_safety's sandbox and
#: deny-list rules, which name those two only.
_CHECK_DISABLE_RE = re.compile(
    r"(?i)\b(?:disable|deactivate|turn\s+off|switch\s+off|skip|bypass"
    r"|circumvent|weaken|remove|silence|suppress|get\s+around|work\s+around"
    r"|abschalten|deaktivier\w*|umgeh\w*|ausschalten)\b"
    r"(?:\s+\w+){0,3}?\s+(?:the\s+|a\s+|any\s+|all\s+|die\s+|den\s+|das\s+)?"
    r"(?:approval|confirm\w*|permission|security|safety|review|lint|test"
    r"|ci|pre-commit|hook|guard|gate|check|codeowners|freigabe\w*"
    r"|pr(?:ü|ue)fung\w*)\w*")

#: The protected zone, by name: the files api_client guards with a confirm
#: (``_DEFAULT_PATH_PROTECTED_GLOBS``) plus the review configuration.
_PROTECTED_NAMES = (
    "api_client.py", "kit_confirm.py", "tab_agent.py", "sandbox.py",
    "kit_settings.py", "user_settings.py", "approval_guard.py",
    "file_confirm.py", "terminal_confirm.py", "approval_answers.py",
    "cli_gate.py", "approval_code.py", "managed_policy.py",
    "principles_addendum.md", "settings.json", "settings.local.json",
    "mcp_servers.json", "trusted_workspaces.json", ".delfin_settings.json",
    "CODEOWNERS", ".github/", ".delfin/skills", "knowledge_check.py",
    "knowledge_share.py",
)
_CHANGE_VERB_RE = re.compile(
    r"(?i)\b(?:edit|modify|change|patch|rewrite|relax|loosen|update|replace"
    r"|delete|remove|overwrite|write\s+to|append\s+to|add\s+to|truncate"
    r"|chmod|ändern|bearbeit\w*|überschreib\w*)\b")

# PUBLIC only ----------------------------------------------------------------
_FILE_EXT = frozenset({
    "py", "md", "txt", "json", "yaml", "yml", "toml", "cfg", "ini", "sh",
    "inp", "out", "xyz", "log", "csv", "html", "js", "ts", "c", "h", "cpp",
    "f90", "gbw", "hess", "pdf", "png", "svg", "ipynb", "lock", "gz", "zip",
})
_FQDN_RE = re.compile(
    r"(?<![\w@/.-])(?:[a-z0-9][a-z0-9-]*\.){2,}[a-z][a-z0-9-]*\b", re.I)
_NODE_RE = re.compile(
    r"\b(?:[a-z]{2,4}\d[a-z]\d{2,4}|[a-z]{2,5}\d{3,5})\b")
_HANDLE_RE = re.compile(r"(?<![\w.@])@[A-Za-z0-9](?:[A-Za-z0-9-]{1,38})\b")
_UNITS = (r"kcal/mol|kj/mol|kJ/mol|kcal|eV|meV|Eh|hartree|Hartree|cm-1"
          r"|cm\^-1|cm⁻¹|nm|Å|A|K|ppm|Hz|MHz|GHz|mV|V|mol/L|mM|µM|%")
_NUMBER_RE = re.compile(
    rf"[-+−]?\d+(?:[.,]\d+)?\s*(?:{_UNITS})(?![\w/])"
    r"|[-+−]?\b\d+\.\d+\b")
_SAMPLE_RE = re.compile(r"\b[A-Z]{1,4}[-_]?\d{2,}[A-Za-z0-9-]*\b")


def _identity_default() -> dict:
    try:
        from .api_client import _identity_for_commit_texts
        return dict(_identity_for_commit_texts())
    except Exception:
        return {}


def _find(rule: str, reason: str, pattern: re.Pattern, text: str,
          out: list, *, skip=None, first_only: bool = False) -> None:
    for m in pattern.finditer(text):
        snippet = m.group(0).strip()
        if not snippet or (skip and skip(snippet, m)):
            continue
        start = m.start() + (len(m.group(0)) - len(m.group(0).lstrip()))
        out.append(Finding(rule, reason, snippet, start, start + len(snippet)))
        if first_only:
            return


def _span_of(text: str, snippet: str) -> tuple[int, int]:
    i = text.find(snippet)
    return (i, i + len(snippet)) if i >= 0 else (-1, -1)


def _secret_findings(text: str, out: list) -> None:
    """Where DELFIN's own redaction would change the text."""
    try:
        from .output_guard import scrub_secrets
        scrubbed = scrub_secrets(text)
    except Exception as exc:
        out.append(Finding("secret", f"secret scan unavailable ({exc!r}); "
                           "nothing is offered unscanned", "", -1, -1))
        return
    if scrubbed == text:
        return
    sm = difflib.SequenceMatcher(a=text, b=scrubbed, autojunk=False)
    for tag, a0, a1, _b0, _b1 in sm.get_opcodes():
        if tag in ("replace", "delete") and a1 > a0:
            out.append(Finding(
                "secret", "a credential or token DELFIN's redaction removes",
                text[a0:a1][:80], a0, a1))


def _instruction(text: str, out: list) -> None:
    """Instructions about permissions, checks or protected files."""
    try:
        from .skill_safety import _TEXT_RULES, _folded
    except Exception as exc:
        out.append(Finding("instruction", f"instruction scan unavailable "
                           f"({exc!r})", "", -1, -1))
        return
    folded = _folded(text)
    for name, pattern, _template in _TEXT_RULES:
        # A URL is judged by the URL rules below: a citation is data in a
        # lesson, while skill_safety reads it as a fetch a recipe performs.
        if name in ("network-detour", "deletion", "scheduler-override"):
            continue
        m = pattern.search(text) or (pattern.search(folded)
                                     if folded != text else None)
        if m:
            snippet = m.group(0).strip()[:80]
            s, e = _span_of(text, snippet)
            out.append(Finding(name, "an instruction about permissions, "
                               "checks or the sandbox", snippet, s, e))
    m = _CHECK_DISABLE_RE.search(text) or (
        _CHECK_DISABLE_RE.search(folded) if folded != text else None)
    if m:
        snippet = m.group(0).strip()[:80]
        s, e = _span_of(text, snippet)
        out.append(Finding("check-disable", "tells the reader to switch off "
                           "a check, a gate or a review", snippet, s, e))
    try:
        from .untrusted import flags
        for flag in flags(folded):
            kind, _, phrase = flag.partition(": ")
            s, e = _span_of(text, phrase)
            out.append(Finding("instruction", f"instruction-shaped text "
                               f"({kind})", phrase, s, e))
    except Exception:
        pass
    for sentence in re.split(r"(?<=[.!?;])\s+|\n", text):
        names = [n for n in _PROTECTED_NAMES if n in sentence]
        if names and _CHANGE_VERB_RE.search(sentence):
            s, e = _span_of(text, sentence.strip())
            out.append(Finding("protected-zone", "names a protected file "
                               "together with a change to it",
                               sentence.strip()[:80], s, e))


def instruction_findings(text: str) -> list[Finding]:
    """The instruction subset of the check, without personal-data rules.

    Used where a text is refused for what it TELLS the reader, not for
    what it reveals: a swarm entry at write time and every repository
    entry at load time (defence in depth behind the share path)."""
    out: list[Finding] = []
    _instruction(str(text or ""), out)
    return _dedupe(out)


def _dedupe(findings: list[Finding]) -> list[Finding]:
    seen: set = set()
    out = []
    for f in findings:
        key = (f.rule, f.snippet, f.start)
        if key not in seen:
            seen.add(key)
            out.append(f)
    return out


def check(text: str, *, public: bool, identity: Optional[dict] = None,
          names: Iterable[str] = (), releasable: Iterable[str] = (),
          unpublished_terms: Optional[Iterable[str]] = None) -> list[Finding]:
    """Every reason not to offer ``text``; empty when it may be offered.

    ``unpublished_terms`` defaults to ``agent.knowledge.unpublished_terms``
    from the user's settings. See the module docstring for the rules."""
    text = str(text or "")
    who = _identity_default() if identity is None else dict(identity)
    out: list[Finding] = []

    try:
        from .skill_share import _control_char_findings, _machine_ref_snippets
        for msg in _control_char_findings(text, "lesson"):
            out.append(Finding("control-character", msg, "", -1, -1))
        for _line_no, snippet, kind in _machine_ref_snippets(text):
            s, e = _span_of(text, snippet)
            rule = "email" if kind == "user@host" and _EMAIL_RE.fullmatch(
                snippet) else "machine-reference"
            out.append(Finding(rule, f"{kind} names a machine or a person",
                               snippet, s, e))
    except Exception as exc:
        out.append(Finding("machine-reference", f"scan unavailable "
                           f"({exc!r})", "", -1, -1))
    _secret_findings(text, out)
    _find("email", "an email address is personal data", _EMAIL_RE, text, out)
    try:
        from .output_guard import home_paths_in
        for hit in home_paths_in(text, home=who.get("home", ""),
                                 account=who.get("account", ""),
                                 host=who.get("host", "")):
            s, e = _span_of(text, hit)
            out.append(Finding("personal-path", "a home path, account or "
                               "machine name of this user", hit, s, e))
    except Exception:
        pass
    _find("licensed-marker", "marks licensed, internal or unpublished "
          "material", _LICENSED_RE, text, out)
    _find("refcode-like", "has the shape of a CSD refcode", _REFCODE_RE,
          text, out, skip=lambda s, _m: s.rstrip("0123456789")
          in _REFCODE_WORDS)
    terms = unpublished_terms
    if terms is None:
        try:
            from delfin.user_settings import load_settings
            cfg = ((load_settings() or {}).get("agent") or {}).get(
                "knowledge") or {}
            terms = cfg.get("unpublished_terms") or ()
        except Exception:
            terms = ()
    for term in terms or ():
        term = str(term or "").strip()
        if len(term) >= 3:
            _find("unpublished-term", "a term you flagged as unpublished",
                  re.compile(re.escape(term), re.I), text, out)
    _find("user-specific", "written about one person or session; write the "
          "general lesson instead", _USER_SPECIFIC_RE, text, out,
          first_only=True)
    for m in _URL_RE.finditer(text):
        url = m.group(0).rstrip(".,;:")
        host = re.sub(r"^[a-z][a-z0-9+.-]*://", "", url, flags=re.I)
        host = host.split("/", 1)[0].split("@")[-1]
        port = ""
        if ":" in host:
            host, port = host.rsplit(":", 1)
        if _INTERNAL_HOST_RE.match(host) or (port and port not in
                                             ("80", "443")):
            out.append(Finding("internal-url", "an internal address",
                               url, m.start(), m.start() + len(url)))
        elif public and not any(host == h or host.endswith("." + h)
                                for h in _PUBLIC_URL_HOSTS):
            out.append(Finding("url", "a URL outside the public allow-list "
                               f"({', '.join(_PUBLIC_URL_HOSTS)})",
                               url, m.start(), m.start() + len(url)))
    _instruction(text, out)
    # The peer-decision shapes (#269): a shared lesson that reads like an
    # approval or a verdict would reach every reader as one.
    try:
        from .untrusted import _DECISION_RE
        m = _DECISION_RE.search(text)
        if m:
            out.append(Finding("decision-shape", "reads like an approval or "
                               "a verdict; a lesson states what to do and "
                               "why, never that something is approved",
                               m.group(0), m.start(), m.end()))
    except Exception:
        pass

    if public:
        url_spans = [(m.start(), m.end()) for m in _URL_RE.finditer(text)]

        def _in_url(m) -> bool:
            return any(a <= m.start() < b for a, b in url_spans)

        def _not_host(snippet, m) -> bool:
            last = snippet.rsplit(".", 1)[-1].lower()
            return _in_url(m) or last in _FILE_EXT or bool(
                _EMAIL_RE.search(text[max(0, m.start() - 64):m.end()]))
        _find("hostname", "a host name", _FQDN_RE, text, out, skip=_not_host)
        _find("hostname", "looks like a cluster node name", _NODE_RE, text,
              out, skip=lambda _s, m: _in_url(m))
        _find("username", "an @handle names a person", _HANDLE_RE, text, out,
              skip=lambda _s, m: _in_url(m))
        _find("number", "a value, possibly from an internal calculation; "
              "mark it releasable if it is published", _NUMBER_RE, text, out,
              skip=lambda _s, m: _in_url(m))
        _find("sample-id", "looks like a sample or project identifier",
              _SAMPLE_RE, text, out, skip=lambda _s, m: _in_url(m))
        for name in names or ():
            name = str(name or "").strip()
            if len(name) >= 3:
                _find("project-name", "the name of the project the lesson "
                      "came from", re.compile(
                          rf"(?<![\w-]){re.escape(name)}(?![\w-])", re.I),
                      text, out)

    released = {str(r).strip() for r in releasable or () if str(r).strip()}
    out = [f for f in out if not (f.releasable and f.snippet in released)]
    return _dedupe(out)


def render(findings: Iterable[Finding]) -> str:
    """One line per finding, for the dialog and the CLI."""
    lines = [f"  - {f.line()}" for f in findings]
    return "\n".join(lines) if lines else "  (no findings)"
