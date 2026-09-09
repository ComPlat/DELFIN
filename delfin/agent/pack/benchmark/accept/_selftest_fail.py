"""Acceptance self-test: always fails, with a reason on stdout."""
import sys
print("this is the reason the artifact was rejected")
sys.exit(1)
