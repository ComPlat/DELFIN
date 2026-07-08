from pathlib import Path

import pytest

from delfin.dashboard.tab_literature import (
    _resolve_literature_path,
    _safe_literature_entry_name,
)


@pytest.mark.parametrize(
    "name",
    [
        "../evil.txt",
        "subdir/paper.pdf",
        r"subdir\paper.pdf",
        "/tmp/evil.txt",
        "C:/tmp/evil.txt",
        ".",
        "..",
        "bad\nname.txt",
        "bad\x00name.txt",
    ],
)
def test_literature_entry_name_rejects_paths_and_control_chars(name):
    with pytest.raises(ValueError):
        _safe_literature_entry_name(name)


def test_literature_entry_name_accepts_plain_names():
    assert _safe_literature_entry_name(" Paper 1.pdf ") == "Paper 1.pdf"


def test_literature_path_must_remain_inside_root(tmp_path):
    root = tmp_path / "literature"
    root.mkdir()
    inside = root / "paper.pdf"
    inside.write_text("ok", encoding="utf-8")

    assert _resolve_literature_path(root, "paper.pdf") == inside.resolve()

    with pytest.raises(ValueError):
        _resolve_literature_path(root, "../outside.txt")


def test_literature_path_blocks_symlink_escape(tmp_path):
    root = tmp_path / "literature"
    outside = tmp_path / "outside"
    root.mkdir()
    outside.mkdir()
    link = root / "link"
    try:
        link.symlink_to(outside, target_is_directory=True)
    except OSError:
        pytest.skip("symlinks unavailable")

    with pytest.raises(ValueError):
        _resolve_literature_path(root, Path("link") / "payload.txt")
