#!/usr/bin/env python3
"""Wrapper: chem_opt_is_a_minimum.py for chloronitrobenzene."""
import importlib.util
import sys
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "_chem_acc", Path(__file__).parent / "chem_opt_is_a_minimum.py")
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
sys.exit(_mod.main([__file__, sys.argv[1] if len(sys.argv) > 1 else "", "chloronitrobenzene"]))
