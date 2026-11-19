"""Verify the website pages: anchors, assets and pinned source links.

Run from the repository root: `python3 -P website/check_pages.py`.
No dependencies, no network, no browser. Checks, in order of what would
embarrass us most:

1. every pinned GitHub source link points at a path that EXISTS at the
   revision the link names (the page's whole claim to being source-linked);
2. every internal #anchor resolves to an id on the page it names;
3. every local asset referenced exists;
4. the stylesheets and the script are syntactically plausible.
"""
import html.parser
import re
import subprocess
import sys
from pathlib import Path

WEB = Path("website")
PAGES = ["framework.html", "index.html"]


class Collect(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids: set[str] = set()
        self.links: list[str] = []
        self.assets: list[str] = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if a.get("id"):
            self.ids.add(a["id"])
        for key in ("href", "src"):
            value = a.get(key) or ""
            if not value:
                continue
            if value.startswith(("http://", "https://", "mailto:", "data:")):
                self.links.append(value)
            elif value.startswith("#") or ".html#" in value:
                self.links.append(value)
            else:
                self.assets.append(value.split("#", 1)[0])


def main() -> int:
    bad: list[str] = []
    parsed: dict[str, Collect] = {}
    for name in PAGES:
        path = WEB / name
        if not path.is_file():
            bad.append(f"missing page: {name}")
            continue
        c = Collect()
        c.feed(path.read_text(encoding="utf-8"))
        parsed[name] = c

    # 1. pinned source links
    pin_re = re.compile(
        r"https://github\.com/ComPlat/DELFIN/(?:blob|tree)/([0-9a-f]{7,40})/([^\"#?]+)")
    pins = 0
    for name, c in parsed.items():
        for link in c.links:
            m = pin_re.match(link)
            if not m:
                continue
            rev, target = m.group(1), m.group(2).rstrip("/")
            pins += 1
            probe = subprocess.run(
                ["git", "cat-file", "-e", f"{rev}:{target}"],
                capture_output=True)
            if probe.returncode != 0:
                bad.append(f"{name}: {target} does not exist at {rev[:9]}")
    print(f"pinned source targets checked: {pins}")

    # 2. internal anchors
    for name, c in parsed.items():
        for link in c.links:
            if link.startswith("#"):
                if link[1:] and link[1:] not in c.ids:
                    bad.append(f"{name}: no id for {link}")
            elif ".html#" in link:
                page, _, frag = link.partition("#")
                target = parsed.get(page)
                if target is None:
                    bad.append(f"{name}: links into unparsed {page}")
                elif frag and frag not in target.ids:
                    bad.append(f"{name}: {page} has no id {frag}")

    # 3. local assets
    for name, c in parsed.items():
        for asset in c.assets:
            if not (WEB / asset).exists():
                bad.append(f"{name}: missing asset {asset}")

    # 4. stylesheets and script, plausibly
    for css in ("styles.css", "framework.css"):
        text = (WEB / css).read_text(encoding="utf-8")
        if text.count("{") != text.count("}"):
            bad.append(f"{css}: unbalanced braces")
    js = (WEB / "framework.js").read_text(encoding="utf-8")
    for opener, closer in (("{", "}"), ("(", ")"), ("[", "]")):
        if js.count(opener) != js.count(closer):
            bad.append(f"framework.js: unbalanced {opener}{closer}")

    if bad:
        print("\nFAILED:")
        for line in bad:
            print("  " + line)
        return 1
    print("anchors, assets, stylesheets and script: all consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
