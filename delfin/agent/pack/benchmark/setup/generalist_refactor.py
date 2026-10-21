#!/usr/bin/env python3
"""Setup for generalist_refactor: rename + move across a package.

``shoplib/`` -- ``calc_total`` lives in shoplib/utils.py beside two
helpers that stay there (``format_price``, ``clamp``). It is reached
three ways: ``from .utils import calc_total`` (cart, checkout),
``utils.calc_total`` through a module import (reports), and the package
export in ``__init__``; tests import it both ways, and docs/usage.md
names it. 19 occurrences in 8 files, test names included.
Precondition: the suite is green and the reference count is as stated.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (GITIGNORE_PY, commit, expect, main, mark_base,
                             new_repo, pytest_rc, write)

FILES = {
    ".gitignore": GITIGNORE_PY,
    "README.md": "# shoplib\n\nA small shop backend. See docs/usage.md.\n",
    "pyproject.toml": ("[project]\nname = \"shoplib\"\nversion = \"2.0.0\"\n\n"
                       "[tool.pytest.ini_options]\npythonpath = [\".\"]\n"
                       "testpaths = [\"tests\"]\n"),
    "docs/usage.md": (
        "# Usage\n\n```python\nfrom shoplib import calc_total\n\n"
        "calc_total([(250, 2), (100, 1)], tax_rate=0.19)  # 714\n```\n\n"
        "`calc_total(items, tax_rate=0.0)` takes (price_cents, quantity)\n"
        "pairs and returns the gross total in cents, rounded half up.\n"),
    "shoplib/__init__.py": (
        "\"\"\"A small shop backend.\"\"\"\n\n"
        "from .utils import calc_total, format_price\n\n"
        "__all__ = [\"calc_total\", \"format_price\"]\n"),
    "shoplib/utils.py": (
        "\"\"\"Shared helpers.\"\"\"\n\n"
        "from decimal import ROUND_HALF_UP, Decimal\n\n\n"
        "def calc_total(items, tax_rate=0.0):\n"
        "    \"\"\"Gross total in cents of (price_cents, quantity) pairs.\"\"\"\n"
        "    net = sum(price * qty for price, qty in items)\n"
        "    gross = Decimal(net) * (Decimal(1) + Decimal(str(tax_rate)))\n"
        "    return int(gross.quantize(Decimal(1), rounding=ROUND_HALF_UP))\n\n\n"
        "def format_price(cents):\n"
        "    \"\"\"12345 -> '123.45 EUR'.\"\"\"\n"
        "    return f\"{cents // 100}.{cents % 100:02d} EUR\"\n\n\n"
        "def clamp(value, low, high):\n"
        "    \"\"\"``value`` limited to [low, high].\"\"\"\n"
        "    return max(low, min(high, value))\n"),
    "shoplib/cart.py": (
        "\"\"\"A shopping cart.\"\"\"\n\n"
        "from .utils import calc_total, clamp\n\n\n"
        "class Cart:\n"
        "    def __init__(self):\n"
        "        self.items = []\n\n"
        "    def add(self, price_cents, qty=1):\n"
        "        self.items.append((price_cents, clamp(qty, 1, 99)))\n\n"
        "    def total(self, tax_rate=0.0):\n"
        "        return calc_total(self.items, tax_rate)\n"),
    "shoplib/checkout.py": (
        "\"\"\"Checkout.\"\"\"\n\n"
        "from .utils import calc_total, format_price\n\n\n"
        "def receipt(items, tax_rate):\n"
        "    \"\"\"One line with the gross total.\"\"\"\n"
        "    total = calc_total(items, tax_rate)\n"
        "    return f\"Total: {format_price(total)}\"\n"),
    "shoplib/reports.py": (
        "\"\"\"Sales reports.\"\"\"\n\n"
        "from . import utils\n\n\n"
        "def daily_revenue(orders, tax_rate=0.0):\n"
        "    \"\"\"Sum of the gross totals of ``orders`` (lists of items).\"\"\"\n"
        "    return sum(utils.calc_total(o, tax_rate) for o in orders)\n"),
    "shoplib/discounts.py": (
        "\"\"\"Discounts.\"\"\"\n\n"
        "from .utils import clamp\n\n\n"
        "def percent_off(cents, percent):\n"
        "    \"\"\"``cents`` reduced by ``percent`` (0-100), rounded down.\"\"\"\n"
        "    return cents * (100 - clamp(percent, 0, 100)) // 100\n"),
    "tests/test_utils.py": (
        "from shoplib.utils import calc_total, clamp, format_price\n\n\n"
        "def test_calc_total_without_tax():\n"
        "    assert calc_total([(250, 2), (100, 1)]) == 600\n\n\n"
        "def test_calc_total_rounds_half_up():\n"
        "    assert calc_total([(250, 2), (100, 1)], tax_rate=0.19) == 714\n\n\n"
        "def test_format_price():\n"
        "    assert format_price(12345) == \"123.45 EUR\"\n\n\n"
        "def test_clamp():\n"
        "    assert clamp(150, 1, 99) == 99\n"),
    "tests/test_cart.py": (
        "from shoplib import calc_total\nfrom shoplib.cart import Cart\n\n\n"
        "def test_cart_total_matches_calc_total():\n"
        "    c = Cart()\n    c.add(250, 2)\n    c.add(100)\n"
        "    assert c.total() == calc_total([(250, 2), (100, 1)])\n"),
    "tests/test_checkout.py": (
        "from shoplib.checkout import receipt\n\n\n"
        "def test_receipt():\n"
        "    assert receipt([(1000, 1)], 0.19) == \"Total: 11.90 EUR\"\n"),
    "tests/test_reports.py": (
        "from shoplib.reports import daily_revenue\n\n\n"
        "def test_daily_revenue():\n"
        "    assert daily_revenue([[(100, 1)], [(200, 2)]]) == 500\n"),
    "tests/test_discounts.py": (
        "from shoplib.discounts import percent_off\n\n\n"
        "def test_percent_off():\n"
        "    assert percent_off(1000, 15) == 850\n"),
}


def build(workspace: str):
    repo = new_repo(workspace, "shoplib")
    write(repo, FILES)
    mark_base(repo, commit(repo, "shoplib 2.0"))
    rc, out = pytest_rc(repo)
    expect("the given suite is green", rc == 0, out)
    refs = sum(text.count("calc_total") for text in FILES.values())
    files = sum(1 for text in FILES.values() if "calc_total" in text)
    expect("19 occurrences in 8 files", (refs, files) == (19, 8),
           f"{refs} references in {files} files")
    return repo


if __name__ == "__main__":
    main(build)
