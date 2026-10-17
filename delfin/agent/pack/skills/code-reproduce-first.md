---
name: code-reproduce-first
description: A failing run or test BEFORE any change, in a repository the user pointed at
domains: code
---
# Reproduce before you fix

In a foreign repository you know nothing about. The report may name the
symptom wrongly; only a failing thing you have run yourself is evidence.

## Steps
1. Read the report and the project layout (README, `ls`, entry point).
2. Run the offered command or test and keep its output. If nothing
   fails, write or run the test that would catch the reported
   behaviour and show it fail.
3. Say in one line what fails, with the actual output.

## If you cannot reproduce
Stop and say exactly what you ran and what you saw. Do not patch on a
guess; ask what differs on the reporter's side.

## Never
- Edit code before a failing run or test exists.
- Treat a passing suite as a reproduction of nothing — if the bug does
  not show, the reproduction is missing, not done.
