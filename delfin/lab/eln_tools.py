"""The ELN tools: four reads, two draft writes a person approves.

Input: a tool name, its arguments as a model sent them, and the session's
``KitToolPermissions``. Output: a JSON string, like every other lab tool
(:mod:`delfin.lab.tools`).

==========  ==================================================================
eln_find    one collection, or the samples / reactions in it, or the ones a
            question selects -- a READ, nothing changes on the server
eln_read    one sample or reaction by id with its analyses; a READ
eln_pending the writes waiting for a person's approval, in words they can
            check; a READ of the lab's own ledger
eln_draft   creates ONE element as a DRAFT; filed for approval, and only a
            person's approval makes the second call send
eln_attach  attaches ONE file (a run's export, a spectrum) to ONE element,
            in the exporter's chunked sequence; filed for approval the same
            way, and the grant covers the whole file by its sha256
eln_approve REFUSES an agent session. The person's command is
            ``delfin-agent eln approve <id>``; this name exists so a model
            that tries it is told which route is the person's, rather than
            accidentally reaching a helper that would refuse it anyway.
==========  ==================================================================

WHERE THE PERMISSION MODE IS AND IS NOT READ. It is read at exactly one
point -- the refusal in :func:`execute` for a write that has no approval --
and only to SAY that it decides nothing. It never reaches
:mod:`delfin.lab.eln_write`, which takes no mode parameter, so a call in
``bypassPermissions`` files the same approval request as a call in ``plan``.

WHAT A TOOL CANNOT DO. It cannot choose the server: the base URL is the
credential's companion setting (:func:`base_url`), not an argument, so a
prompt that carries a URL cannot redirect a write. It cannot pass a
username, a password or a token. It cannot list an element type the route
table does not cover (``eln_routes`` refuses the (verb, path), and the
refusal is quoted back).
"""

from __future__ import annotations

import fnmatch
import json
import os
from typing import Optional

from . import eln_api, eln_write

#: The names this module answers. Keep in step with
#: ``api_client._LAB_TOOL_NAMES`` and the catalogue schemas; the test
#: ``tests/test_eln_tools_are_wired.py`` holds the three together.
TOOL_NAMES = frozenset({"eln_find", "eln_read", "eln_pending", "eln_draft",
                        "eln_approve", "eln_attach"})

#: Read names only. These need no approval and never write.
READ_NAMES = frozenset({"eln_find", "eln_read", "eln_pending"})

#: The element types a read may ask for. The plural the route spells.
READ_ELEMENTS = ("samples", "reactions", "research_plans", "wellplates",
                 "screens")

#: The settings key the base URL is read from, and the environment variable
#: that overrides it. A URL, never a secret.
BASE_URL_SETTING = "lab.eln.base_url"
BASE_URL_ENV = "DELFIN_ELN_URL"


def base_url() -> str:
    """The ELN's address for this installation, or "".

    Read from the environment first (a per-process override for a fake
    server or a throwaway Docker instance), then the user's settings. The
    value is a URL: the credential travels in the Authorization header,
    never in the address, so this is safe to show.
    """
    env = os.environ.get(BASE_URL_ENV, "").strip()
    if env:
        return env.rstrip("/")
    try:
        from delfin.user_settings import load_settings
        settings = load_settings() or {}
        got = ((settings.get("lab") or {}).get("eln") or {}).get("base_url")
        return str(got or "").strip().rstrip("/")
    except Exception:                                   # noqa: BLE001
        return ""


def _credential(name: str) -> str:
    """The credential store's answer, by NAME. The only route to a secret."""
    from delfin.agent.credentials import load_credential
    return load_credential(name) or ""


def _client(transport=None) -> eln_api.Client:
    url = base_url()
    if not url:
        raise eln_api.ClientError(
            f"no ELN address configured; set {BASE_URL_SETTING} (environment "
            f"{BASE_URL_ENV} overrides it). The credential is separate: store "
            f"it as {eln_api.TOKEN_NAME}.")
    return eln_api.Client(url, _credential, transport=transport)


def _dump(obj) -> str:
    return json.dumps(obj, default=str, indent=1, sort_keys=True)


def _refuse_secret_path(path: str) -> Optional[dict]:
    """The agent's own deny list, applied before this tool reads anything.

    ``eln_attach`` reads a path off the tool call. The agent's file tools
    already refuse a deny-listed path (``*.key``, ``.ssh/**``,
    ``credentials*``, ``.netrc``, ``.env*`` ...), and that gate is the whole
    reason a key cannot be pulled into a conversation. A second reader without
    it is that gate missing, so this asks the same list rather than growing a
    second one: the globs come from ``api_client``, and the check is done on
    the path AS GIVEN and on its resolved form, because ``~/.ssh/id_rsa`` and
    its target are the same file under a different spelling.

    Returns the refusal dict, or None when the path is not deny-listed. A path
    that cannot be resolved is refused too: "I could not tell" is not a pass.
    """
    from pathlib import Path as _Path

    from delfin.agent.api_client import _DEFAULT_PATH_DENY_GLOBS as globs

    raw = str(path)
    try:
        resolved = _Path(raw).expanduser().resolve()
    except Exception:                                   # noqa: BLE001
        return {"error": f"cannot resolve {raw!r}, so it is not read",
                "code": "refused_path"}
    for candidate in (raw, str(resolved), resolved.name):
        rel = candidate.replace("\\", "/")
        low = rel.lower()
        for g in globs:
            if fnmatch.fnmatch(rel, g) or fnmatch.fnmatchcase(low, g.lower()):
                return {
                    "error": (f"{raw!r} matches the deny list entry {g!r}; "
                              "this tool does not read it. Nothing was read "
                              "and nothing was filed. Attach a file the lab "
                              "actually produced (a run export, a spectrum) "
                              "rather than a credential or a key."),
                    "code": "refused_path",
                    "matched": g,
                }
    return None


# ---------------------------------------------------------------------------
# The questions
# ---------------------------------------------------------------------------

def _find(args: dict, transport=None) -> dict:
    """A collection with its contents, or a search across one.

    Two shapes, because a chemist asks two things: "what is in collection
    X" (``collection``) and "which samples in X have no analysis" (a
    filter, encoded as the search body the server takes).
    """
    client = _client(transport)
    collection = str(args.get("collection") or "").strip()
    element = str(args.get("element") or "").strip() or "samples"
    query = args.get("query")
    if isinstance(query, str) and query.strip():
        try:
            query = json.loads(query)
        except ValueError:
            return {"error": f"query must be a JSON object, not {query!r}"}
    if query is not None and not isinstance(query, dict):
        return {"error": "query must be a JSON object"}

    out: dict = {"element": element}
    if collection:
        row = eln_api.resolve_collection(client, collection)
        out["collection"] = {"id": row.get("id"),
                             "label": row.get("label") or row.get("name")}
        body = dict(query or {})
        body.setdefault("collection_id", row.get("id"))
    else:
        body = dict(query or {})
    rows = client.search_all(element, body, max_pages=int(
        args.get("max_pages") or 20))
    out["rows"] = _condense(element, rows)
    out["count"] = len(out["rows"])
    out["population"] = (
        f"{len(out['rows'])} {element} the server returned for this search"
        + (f" in collection {out['collection']['label']!r}"
           if out.get("collection") else ""))
    if not rows:
        out["note"] = ("no row matched. The search body the server got was "
                       f"{body!r}; an empty answer can also mean the element "
                       "is not visible to this account")
    return out


def _condense(element: str, rows: list) -> list:
    """The fields a chemist asks about, per row, from what the entity exposes.

    Read from the entities' own ``expose`` lists
    (:data:`eln_api.SOURCES`): a sample carries ``name``, ``short_label``,
    ``external_label``, ``sum_formula`` and a ``container`` whose children
    are its analyses; a reaction carries ``name``, ``short_label``,
    ``status`` and its material lists. Everything else is dropped: a
    transcript full of molfiles is a transcript nobody reads.
    """
    out: list = []
    for row in rows:
        if not isinstance(row, dict):
            out.append({"row": str(row)})
            continue
        item = {
            "id": row.get("id"),
            "name": row.get("name") or row.get("short_label"),
            "short_label": row.get("short_label"),
            "external_label": row.get("external_label"),
            "type": row.get("type"),
        }
        if element.startswith("sample"):
            item["sum_formula"] = row.get("sum_formula")
            names = _analysis_names(row)
            # None is a THIRD answer, not an absence: the row does not say
            # whether this sample has an analysis, and reporting it as
            # "has none" would answer the chemist's question wrongly.
            item["analyses"] = names
            item["has_analysis"] = None if names is None else bool(names)
        elif element.startswith("reaction"):
            item["status"] = row.get("status")
            for key in ("reactants", "products", "solvents",
                        "starting_materials"):
                value = row.get(key)
                if isinstance(value, list):
                    item[key] = [_material_name(m) for m in value]
        else:
            item["body"] = bool(row.get("body"))
        out.append(item)
    return out


def _material_name(material) -> str:
    if not isinstance(material, dict):
        return str(material)
    for key in ("name", "short_label", "id"):
        if material.get(key) is not None:
            return str(material[key])
    return "?"


def _analysis_names(row: dict) -> Optional[list]:
    """The analysis kinds a sample carries, or None when the row cannot say.

    ``Entities::SampleEntity`` exposes ``container`` at level 0, so a row
    that carries no ``container`` key was serialised for a LIST
    (``displayed_in_list``) or withheld it -- and an empty list then means
    "no analyses", which is the one answer a chemist must not be given
    wrongly: "which samples have no analysis" would name every sample in a
    list answer.

    So there are three outcomes, and they are not the same answer:
    ``None`` = this row does not say (read the element itself); ``[]`` = it
    says there are none; a list of kinds = these are its analyses.

    An analysis is a child container with ``container_type == 'analysis'``
    whose ``extended_metadata.kind`` names the technique
    (``container_entity.rb``'s ``expose :children`` and its
    ``extended_metadata[:kind]`` mapping).
    """
    if "container" not in row:
        return None
    container = row.get("container")
    if not isinstance(container, dict):
        return None
    names: list = []
    for child in container.get("children") or []:
        if not isinstance(child, dict):
            continue
        kind = ""
        meta = child.get("extended_metadata")
        if isinstance(meta, dict):
            kind = str(meta.get("kind") or "")
        if child.get("container_type") == "analysis" or kind:
            names.append(kind or str(child.get("name") or "analysis"))
    return names


def _read(args: dict, transport=None) -> dict:
    """One element by id, with its analyses spelled out."""
    client = _client(transport)
    element = str(args.get("element") or "").strip() or "samples"
    element_id = args.get("id")
    if element_id in (None, ""):
        return {"error": "eln_read needs an id; eln_find lists the ids in a "
                         "collection"}
    row = client.element(element, element_id)
    out = {"element": element, "id": element_id, "row": row}
    if element.startswith("sample"):
        names = _analysis_names(row)
        out["analyses"] = names
        out["has_analysis"] = None if names is None else bool(names)
        if names is None:
            out["note"] = ("this answer carries no container, so it does not "
                           "say whether the sample has an analysis")
    return out


def _pending(args: dict, transport=None) -> dict:
    """The writes waiting for a person, with what each would write."""
    rows = []
    for rec in eln_write.pending():
        rows.append({
            "id": rec.get("id"),
            "what": eln_write.describe(
                rec.get("element", ""), {"name": rec.get("name"),
                                         "collection_id": rec.get(
                                             "collection_id")}),
            "server": rec.get("base_url"),
            "expires_in_s": max(0, int(rec.get("expires", 0)
                                       - eln_write._now())),
            "approve_with": f"delfin-agent eln approve {rec.get('id')}",
        })
    return {"pending": rows, "count": len(rows),
            "note": ("approving is a person's action; an agent session "
                     "cannot do it, and no permission mode can")}


def _draft(args: dict, permissions=None, transport=None) -> dict:
    """Create one DRAFT, or file the approval request that would allow it."""
    client = _client(transport)
    element = str(args.get("element") or "").strip()
    fields = dict(args.get("fields") or {})
    if not isinstance(fields, dict) or not fields:
        # Accept the flat spelling a model tends to send.
        fields = {k: v for k, v in args.items()
                  if k in ("name", "collection_id", "materials", "body",
                           "description", "observation")}
    collection = args.get("collection")
    if collection not in (None, "") and fields.get("collection_id") is None:
        row = eln_api.resolve_collection(client, str(collection))
        fields["collection_id"] = row.get("id")
    mode = str(getattr(permissions, "mode", "") or "")
    session = str(getattr(permissions, "task_session_id", "") or "")
    try:
        return eln_write.create_draft(
            client, element, fields, approval=str(args.get("approval") or ""),
            session=session)
    except eln_write.WriteRefused as exc:
        out = {"error": str(exc), "code": exc.code, "draft": exc.draft,
               "request_id": exc.request_id, "digest": exc.digest}
        if exc.code == "needs_approval":
            # The one place the mode is read, and only to say it decides
            # nothing. A model that reached this in bypassPermissions is
            # told plainly, because "I was allowed to write" is the
            # assumption this tool exists to break.
            out["note"] = (
                f"this is not decided by the permission mode (this call's "
                f"was {mode or 'unset'}), by a prompt, or by a claimed "
                "approval. Only the person's `delfin-agent eln approve "
                f"{exc.request_id}` makes it sendable.")
            out["approve_with"] = f"delfin-agent eln approve {exc.request_id}"
        return out


def _attach(args: dict, transport=None) -> dict:
    """Attach one file to one element, or file the approval that would.

    ``path`` is a file on this machine (a lab run's export, a spectrum);
    ``content`` is accepted as an alternative for a caller that already
    holds the bytes. One of the two. The element type is required because
    the link request names it, and the id because the file has to land
    somewhere.
    """
    path = str(args.get("path") or "").strip()
    content = args.get("content")
    if path:
        # BEFORE the read: this tool is a reader of the machine's files, and
        # the agent's file tools refuse a deny-listed path -- keys, ssh
        # material, credential stores, .env. Hashing the bytes first and filing
        # an approval that shows only a name, a size and a digest would make
        # this route around that gate, not carry it.
        refused = _refuse_secret_path(path)
        if refused is not None:
            return refused
        from pathlib import Path
        try:
            content = Path(path).read_bytes()
        except OSError as exc:
            return {"error": f"cannot read {path!r}: {exc}", "code": "bad_path"}
    if content is None:
        return {"error": "eln_attach needs a path to a file (or content "
                         "bytes); nothing was read and nothing was sent."}
    if isinstance(content, str):
        content = content.encode("utf-8")
    element_type = str(args.get("element_type") or "").strip()
    element_id = args.get("element_id")
    try:
        element_id = int(element_id)
    except (TypeError, ValueError):
        return {"error": "eln_attach needs an integer element_id: the file "
                         "must land on an element, and the id is what names "
                         "it.", "code": "bad_element_id"}
    filename = str(args.get("filename") or "").strip() or (
        path.rsplit("/", 1)[-1] if path else "attachment")
    client = _client(transport)
    try:
        return eln_write.attach_file(
            client, filename=filename, content=content,
            element_type=element_type, element_id=element_id,
            approval=str(args.get("approval") or ""))
    except eln_write.WriteRefused as exc:
        out = {"error": str(exc), "code": exc.code, "draft": exc.draft,
               "request_id": exc.request_id, "digest": exc.digest}
        if exc.code == "needs_approval":
            out["note"] = (
                "this is not decided by the permission mode, by a prompt, or "
                "by a claimed approval. Only the person's `delfin-agent eln "
                f"approve {exc.request_id}` makes it sendable.")
            out["approve_with"] = f"delfin-agent eln approve {exc.request_id}"
        return out


def _approve(args: dict, transport=None) -> dict:
    """Refused: approving is the person's command, never an agent's tool."""
    return {
        "error": ("approving an ELN write is a person's action, through "
                  "`delfin-agent eln approve <id>` in their own terminal; no "
                  "tool call can do it, and this process cannot be made to"),
        "code": "refused_person_only",
        "person_command": "delfin-agent eln approve <id>",
    }


# ---------------------------------------------------------------------------
# The entry point the executor already dispatches to
# ---------------------------------------------------------------------------

def execute(name: str, arguments: dict, permissions=None) -> str:
    """Run one ELN tool. Returns a JSON string, never raises.

    ``transport`` (a test's fake server) is taken from ``arguments`` ONLY
    when it is a callable -- an in-process function, not something a model
    can send over the wire. That is deliberate: it keeps the fake server
    out of the product's surface while letting a test drive the real entry
    point instead of the helper.
    """
    args = arguments if isinstance(arguments, dict) else {}
    transport = args.pop("_transport", None)
    if not callable(transport):
        transport = None
    try:
        if name in ("eln_find", "eln_read", "eln_pending"):
            fn = {"eln_find": _find, "eln_read": _read,
                  "eln_pending": _pending}[name]
            return _dump(fn(args, transport=transport))
        if name == "eln_draft":
            return _dump(_draft(args, permissions=permissions,
                                transport=transport))
        if name == "eln_attach":
            return _dump(_attach(args, transport=transport))
        if name == "eln_approve":
            return _dump(_approve(args))
    except eln_api.ClientError as exc:
        return _dump({"error": str(exc), "status": exc.status,
                      "detail": exc.detail, "route": exc.route})
    except eln_write.WriteRefused as exc:
        return _dump({"error": str(exc), "code": exc.code,
                      "request_id": exc.request_id})
    except Exception as exc:                            # noqa: BLE001
        return _dump({"error": f"{type(exc).__name__}: {exc}"})
    return _dump({"error": f"unknown ELN tool {name!r}"})
