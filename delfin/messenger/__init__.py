"""The messenger gateway: the SimpleX phone side of the agent service.

Layout:
    simplex.py   the daemon connection (local WebSocket API, 127.0.0.1)
    bridge_client.py  the gateway side of the service bridge's JSON lines
    store.py     the pairing state, 0700 under $HOME/.delfin/messenger
    gateway.py   the command loop: only the paired contact, the fixed
                 commands, limits
    __main__.py  pair | run | revoke | status

Everything the phone may say is a fixed set; everything the phone gets
is what the bridge redacted already. Approvals are never granted here:
the gateway names where to approve, the dashboard keeps the right.
"""
