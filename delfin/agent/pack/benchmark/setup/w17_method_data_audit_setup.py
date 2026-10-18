#!/usr/bin/env python3
"""Task setup entry: data audit variant (w17_method_data_audit).

The setup contract (usage / refusing to overwrite / prepared) comes
from generalist_repo.main, like every other benchmark setup.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generalist_repo import main  # noqa: E402
from w17_method_data_repo import build  # noqa: E402

if __name__ == "__main__":
    main(lambda ws: build(ws, variant="audit"))
