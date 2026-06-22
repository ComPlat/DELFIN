# TDD Mode (Test-Driven Development)

Test-first workflow: the Test Agent writes tests based on the acceptance criteria
BEFORE the Builder implements. The Builder then makes the tests pass.

Route: Session Manager → Test Agent (write tests) → Builder (make tests pass) → Reviewer → Test Agent (verify)

Best for:
- Well-defined features with clear acceptance criteria
- Bug fixes where a failing test should be written first
- Refactoring where existing behavior must be preserved
- Any task where "done" can be expressed as passing tests
