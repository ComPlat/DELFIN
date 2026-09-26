"""A refusal is recorded before the client can read the 403.

_refuse recorded after sendall; the serving thread could still be
between the two when the client had already read its answer and looked
at the record (SLURM 7210626: proxy.refused == [] beside a correct 403).
"""

from __future__ import annotations

from delfin.agent import egress_proxy


class _Conn:
    def __init__(self, seen):
        self.seen = seen

    def sendall(self, data):
        # At the moment the answer leaves, the refusal must be on record.
        self.seen.append(("sent", list(self.records)))


def test_the_record_exists_when_the_answer_is_sent():
    records = []
    seen = []
    proxy = egress_proxy.EgressProxy.__new__(egress_proxy.EgressProxy)
    proxy._on_refusal = records.append
    conn = _Conn(seen)
    conn.records = records
    proxy._refuse(conn, 403, "blocked by DELFIN's sandbox: example.org")
    assert seen == [("sent", ["blocked by DELFIN's sandbox: example.org"])]
