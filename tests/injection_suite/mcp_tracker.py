"""A stdio MCP server that plays an issue tracker for one suite case.

Input: newline-delimited JSON-RPC on stdin; argv[1] is a JSON file with the
issues ``search`` returns. Output: answers to ``initialize``,
``tools/list`` (one tool, ``search``, a read-only name) and
``tools/call``. The issue text is the untrusted source of the case; the
server itself does nothing else.
"""

import json
import sys


def main() -> None:
    with open(sys.argv[1], encoding="utf-8") as fh:
        issues = json.load(fh)
    for line in sys.stdin:
        try:
            msg = json.loads(line)
        except ValueError:
            continue
        if "id" not in msg:
            continue
        method = msg.get("method")
        if method == "initialize":
            result = {"protocolVersion": "2025-06-18", "capabilities":
                      {"tools": {}}, "serverInfo": {"name": "tracker",
                                                    "version": "1"}}
        elif method == "tools/list":
            result = {"tools": [{
                "name": "search",
                "description": "Search the project's issues by keyword.",
                "inputSchema": {"type": "object", "properties": {
                    "query": {"type": "string"}}, "required": ["query"]}}]}
        elif method == "tools/call":
            query = str(((msg.get("params") or {}).get("arguments") or {})
                        .get("query") or "").lower()
            hits = [i for i in issues if query in json.dumps(i).lower()]
            result = {"content": [{"type": "text",
                                   "text": json.dumps(hits, indent=1)}]}
        else:
            result = {}
        sys.stdout.write(json.dumps({"jsonrpc": "2.0", "id": msg["id"],
                                     "result": result}) + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
