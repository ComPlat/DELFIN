#!/bin/sh
# The whole fast suite in three passes, run to COMPLETION by one command.
#
# Three passes rather than one because the whole-tree run takes longer than a
# single background slot on this host, and a result nobody finished is not a
# result. Each pass writes its own log under scratch/p3-a1/ and its counts are
# quoted in scratch/PR.md exactly as they appear.
#
# The seven --ignore entries are tests that cannot even be imported here:
# they import `cluster_tunnel`, which is not installed in this container.
# That is a property of the host, not of the change: the same seven fail to
# collect on an untouched origin/main worktree (recorded in scratch/PR.md).
set -u
cd "$(dirname "$0")/../.." || exit 1

IGN="
--ignore=tests/test_a_batch_of_cluster_jobs_is_checked_before_anyone_is_asked.py
--ignore=tests/test_a_batch_waits_for_its_install_without_holding_the_turn.py
--ignore=tests/test_a_cluster_gets_the_pcs_published_commit_installed_beside_the_others.py
--ignore=tests/test_a_cluster_is_added_in_the_window_and_only_proposed_by_an_agent.py
--ignore=tests/test_a_cluster_is_used_only_once_its_host_key_is_pinned.py
--ignore=tests/test_a_cluster_job_comes_back_to_the_session_that_sent_it.py
--ignore=tests/test_the_agent_reaches_a_cluster_only_through_its_tools_and_a_person.py
"

i=1
for chunk in 1 2 3; do
    files=$(ls tests/test_*.py | sed -n "${chunk}~3p")
    echo "=== chunk $chunk: $(echo "$files" | wc -l) files ==="
    PYTHONPATH=. python -m pytest $files $IGN \
        -q -m "not slow" --timeout=120 --timeout-method=signal -p no:randomly \
        > "scratch/p3-a1/suite-third${chunk}.log" 2>&1
    echo "chunk $chunk exit=$?"
    tail -3 "scratch/p3-a1/suite-third${chunk}.log"
    i=$((i + 1))
done
echo "=== done ==="
