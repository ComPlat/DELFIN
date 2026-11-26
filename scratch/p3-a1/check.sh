#!/bin/sh
# The checks this package repeats: the review test, then the suites that
# consume the presets and the slash commands it changes.
set -e
cd "$(dirname "$0")/../.."
PYTHONPATH=. python -m pytest \
  tests/test_a_review_finding_carries_its_location_and_its_failure.py \
  tests/test_sa0_invariants.py tests/test_sa0_preset_fields.py \
  tests/test_subagents.py tests/test_subagent_contracts.py \
  tests/test_subagent_slash_commands.py \
  tests/test_subagent_md_loader.py tests/test_sa0_validate_preset.py \
  tests/test_every_shipped_preset_is_valid.py \
  tests/test_the_routing_table_says_when_not_just_what.py \
  tests/test_a_preset_has_its_tools_rather_than_asking_for_them.py \
  tests/test_slash_commands_introspection.py \
  -q --timeout=180
