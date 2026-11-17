"""The run budget and the delegated-spend ledger, as one object.

Extracted from ``delfin.agent.engine`` (they lived there as ~15 methods)
without any behaviour change: the engine keeps every method name as a
one-line delegate with the same signature, so no caller -- run_limits,
scheduler_daemon, the dashboard or a test -- had to change.

The ledger reads the engine's state (``cost_usd``, the turn counters,
``run_budget_*``, ``token_usage``) through the engine reference it is
given at construction; it does not copy anything. Test engines built with
``AgentEngine.__new__`` never run ``__init__``, so the engine's
delegates install the ledger lazily when it is missing.
"""

from __future__ import annotations

from . import agent_metrics as _agent_metrics

__all__ = ["BudgetLedger"]


class BudgetLedger:
    """The run budget (USD / wall clock / turns), the per-turn cost
    breaker and the delegated-spend metering, reading its state through
    the engine it belongs to.

    Every method below is the code that used to live on AgentEngine,
    moved as-is; ``self._engine`` stands where ``self`` stood.
    """

    def __init__(self, engine) -> None:
        self._engine = engine

    # -- the ceiling dimensions -------------------------------------------

    def _run_budget(self) -> tuple[float, float]:
        """(usd, seconds) run-scoped budget for this SESSION — the cumulative
        ceiling that makes unattended runs deployable (the per-turn cap below
        only stops a single runaway turn; turns compose into unbounded run
        cost without this). 0 disables either dimension. An explicit
        ``self.run_budget_usd``/``self.run_budget_s`` attribute (set e.g. by
        the scheduler daemon per entry) overrides the settings values."""
        engine = self._engine
        usd = float(getattr(engine, "run_budget_usd", 0.0) or 0.0)
        secs = float(getattr(engine, "run_budget_s", 0.0) or 0.0)
        if usd <= 0 or secs <= 0:
            try:
                from delfin.user_settings import load_settings as _ls
                ag = (_ls() or {}).get("agent", {}) or {}
                if usd <= 0:
                    usd = max(0.0, float(ag.get("run_budget_usd", 0.0) or 0.0))
                if secs <= 0:
                    secs = max(0.0, float(ag.get("run_budget_s", 0.0) or 0.0))
            except Exception:
                pass
        return usd, secs

    def _run_ceiling_hit(self, spent_this_turn: float = 0.0) -> str:
        """The run ceiling a RUNNING turn has reached, or "".

        Input: what the running turn has spent so far (the client's
        estimate; ``cost_usd`` holds the finished turns). Output: the
        stop reason the turn ends with -- ``run_budget_s`` once the wall
        clock since ``_run_started_at`` reaches ``run_budget_s``,
        ``run_budget_usd`` once finished plus running spend reaches
        ``run_budget_usd`` -- else "". The turn ceiling is not checked
        here: it counts finished turns, which do not change mid-turn.
        Consulted by the client between tool rounds (``run_ceiling``);
        the 110 % gate before a turn stays. Never raises.
        """
        engine = self._engine
        hit = ""
        try:
            usd, secs = self._run_budget()
            if secs > 0:
                import time as _t
                started = float(getattr(engine, "_run_started_at", 0.0) or 0.0)
                if started > 0 and _t.time() - started >= secs:
                    hit = "run_budget_s"
            if not hit and usd > 0:
                spent = (float(getattr(engine, "cost_usd", 0.0) or 0.0)
                         + max(0.0, float(spent_this_turn or 0.0)))
                if spent >= usd:
                    hit = "run_budget_usd"
        except Exception:
            return ""
        if hit:
            # Once per run (the level dedups): the turn's own notice goes
            # to the stream, this one to whoever is not watching it.
            try:
                from .run_limits import resume_hint
                engine._emit_budget_attention(
                    "ceiling", "Run stopped: budget ceiling reached",
                    f"agent.{hit} was reached inside a turn; the turn "
                    "stopped between two tool rounds, its work is kept and "
                    "the session is saved. To continue: "
                    + resume_hint(engine, hit))
            except Exception:
                pass
        return hit

    def _run_budget_turns(self) -> int:
        """Turn ceiling for this run — 0 disables it.

        The third dimension exists because the first one cannot always be
        measured. A turn is countable whatever the model costs, so this
        is the ceiling that still holds when the USD one is inert (see
        :meth:`_usd_budget_enforced`). Same precedence as the other two:
        an explicit attribute set per scheduler entry beats the setting.
        """
        engine = self._engine
        turns = int(getattr(engine, "run_budget_turns", 0) or 0)
        if turns <= 0:
            try:
                from delfin.user_settings import load_settings as _ls
                ag = (_ls() or {}).get("agent", {}) or {}
                turns = int(ag.get("run_budget_turns", 0) or 0)
            except Exception:
                turns = 0
        return max(0, turns)

    def _turns_taken(self) -> int:
        """Turns this run has actually run, priced or not.

        Deliberately the SUM of the three price states and not the
        unpriced count: a run that starts on an unpriced model and
        switches to a priced one would otherwise stop spending its
        ceiling half way through.
        """
        engine = self._engine
        return (int(getattr(engine, "_measured_cost_turns", 0) or 0)
                + int(getattr(engine, "_non_billing_turns", 0) or 0)
                + int(getattr(engine, "_unpriced_turns", 0) or 0))

    # -- what a turn's cost meant -----------------------------------------

    def _note_turn_price_state(self, cost_delta: float = 0.0) -> None:
        """Count the turn that just ran by what its cost could mean.

        Three ways the number means something, resolved in the order the
        dashboard's own cost line already uses:

        * NON_BILLING — the provider charges no USD for the call at all
          (KIT quota, local hardware). A measured zero.
        * The backend measured it itself: a reported spend above zero, or
          the CLI backend, where the provider states the turn's cost
          (often zero, on a subscription) instead of leaving it to a rate
          table. An observation outranks the table either way.
        * PRICED — a published rate applied to counted tokens.

        Everything else is UNKNOWN: no rate, nothing observed, and a 0.0
        in ``self.cost_usd`` that every USD gate reads as thrift.
        """
        from . import pricing as _pricing
        engine = self._engine
        state = _pricing.resolve(
            str(getattr(engine.client, "model", "") or ""),
            str(getattr(engine, "provider", "") or ""),
        ).state
        if state == _pricing.NON_BILLING:
            engine._non_billing_turns += 1
            engine._last_turn_price_state = "non_billing"
        elif (state == _pricing.PRICED
                or float(cost_delta or 0.0) > 0
                or getattr(engine, "backend", "") == "cli"):
            engine._measured_cost_turns += 1
            engine._last_turn_price_state = "measured"
        else:
            engine._unpriced_turns += 1
            engine._last_turn_price_state = "unknown"

    def _usd_budget_enforced(self) -> bool:
        """Whether the configured USD ceiling actually bounds this run.

        False in one case that matters: a USD budget is set and at least
        one turn ran on a model with no published rate. Such a turn adds
        0.0 to ``self.cost_usd``, so the fraction the gates compare is a
        sum over the turns that happened to be measurable — it can sit at
        0% through a run of any size, and the refusal above 110% can never
        fire. The tier aliases the model picker ships (``opus``,
        ``sonnet``, ``haiku``) are precisely the unpriced ids, so this is
        the default case, not an exotic one.

        A quota-funded or locally served run is NOT this case: its zero is
        a measurement, its ceiling holds, and calling it unmeasured would
        be a false statement in the other direction.
        """
        engine = self._engine
        usd, _secs = self._run_budget()
        return usd > 0 and int(getattr(engine, "_unpriced_turns", 0) or 0) == 0

    def _unmeasured_budget_block(self) -> str:
        """Say ONCE that no monetary ceiling is in force, and what is.

        Silent unless a USD budget is configured AND an unpriced turn has
        run: with no budget there is nothing to mislead anyone about, and
        with every turn priced the ordinary wind-down block does the job.
        Turn counts and token totals are always measured, so they are what
        this offers in place of the ceiling that cannot be enforced.
        """
        engine = self._engine
        usd, secs = self._run_budget()
        unpriced = int(getattr(engine, "_unpriced_turns", 0) or 0)
        if (usd <= 0 or unpriced == 0
                or getattr(engine, "_unmeasured_budget_notice_shown", False)):
            return ""
        engine._unmeasured_budget_notice_shown = True
        turn_cap = self._run_budget_turns()
        # The prompt reaches the model and nobody else. A run budget is
        # what makes an unattended run deployable, so the one fact that
        # it is not in force has to reach whoever is not watching.
        engine._emit_budget_attention(
            "unenforceable",
            "Run budget: the USD ceiling is not in force",
            f"{unpriced} turn(s) ran on a model with no published rate, so "
            f"the ${usd:.2f} ceiling cannot be enforced and neither can the "
            f"per-turn cost breaker."
            + (f" The {turn_cap}-turn ceiling IS in force."
               if turn_cap > 0 else
               " Set agent.run_budget_turns or agent.run_budget_s for a "
               "ceiling that can be."))
        turns = (unpriced + int(getattr(engine, "_measured_cost_turns", 0) or 0)
                 + int(getattr(engine, "_non_billing_turns", 0) or 0))
        tok_in = int(engine.token_usage.get("input", 0) or 0)
        tok_out = int(engine.token_usage.get("output", 0) or 0)
        lines = [
            "# Run budget NOT enforceable in USD (auto-injected, once)",
            f"- {unpriced} of {turns} turn(s) so far ran on a model with no "
            f"published rate, so their cost was never measured. The "
            f"${usd:.2f} ceiling is not in force, and neither is the "
            f"per-turn cost breaker.",
            f"- Measured instead: {turns} turn(s), {tok_in:,} input + "
            f"{tok_out:,} output tokens. Judge the size of this run by "
            f"those, not by a dollar figure.",
        ]
        if turn_cap > 0:
            lines.append(f"- The {turn_cap} turn(s) ceiling IS measured and "
                         "still ends this run.")
        if secs > 0:
            lines.append(f"- The {secs:.0f}s wall-clock ceiling IS measured "
                         "and still ends this run.")
        if turn_cap <= 0 and secs <= 0:
            lines.append("- For a ceiling that can be enforced, set "
                         "agent.run_budget_turns or agent.run_budget_s, or "
                         "pick a model id with a published rate.")
        return "\n".join(lines)

    # -- the fraction spent ------------------------------------------------

    def _run_budget_requests(self) -> int:
        """Request ceiling for this run — 0 disables it. OPTIONAL, default
        OFF: absent from the settings and unset on the engine it adds no
        dimension to :meth:`_run_budget_detail` at all.

        **Cost in prompt text:** when it is in force, the wind-down /
        exhausted budget block and the >110 % refusal name it like the
        other dials — one short setting name in text that is already
        injected every turn past 80 % of any configured budget; with the
        cap off the prompt text is unchanged.

        Same precedence as the other dimensions: an explicit
        ``run_budget_requests`` attribute, which a caller that builds the
        engine itself (a scheduler entry, a test) may set, beats the
        settings value ``agent.run_budget_requests``. Unlike the three
        older dimensions this one has no command-line flag and no
        unattended default in ``run_limits``: it is off until somebody
        asks for it.
        """
        engine = self._engine
        reqs = int(getattr(engine, "run_budget_requests", 0) or 0)
        if reqs <= 0:
            try:
                from delfin.user_settings import load_settings as _ls
                ag = (_ls() or {}).get("agent", {}) or {}
                reqs = int(ag.get("run_budget_requests", 0) or 0)
            except Exception:
                reqs = 0
        return max(0, reqs)

    def _run_budget_tokens(self) -> int:
        """Token ceiling for this run — 0 disables it. OPTIONAL, default
        OFF: absent from the settings and unset on the engine it adds no
        dimension to :meth:`_run_budget_detail` at all.

        **Cost in prompt text:** when it is in force, the wind-down /
        exhausted budget block and the >110 % refusal name it like the
        other dials — one short setting name in text that is already
        injected every turn past 80 % of any configured budget; with the
        cap off the prompt text is unchanged.

        Same precedence as the other dimensions: an explicit
        ``run_budget_tokens`` attribute beats the settings value
        ``agent.run_budget_tokens``. No command-line flag and no
        unattended default: off until somebody asks for it.
        """
        engine = self._engine
        tokens = int(getattr(engine, "run_budget_tokens", 0) or 0)
        if tokens <= 0:
            try:
                from delfin.user_settings import load_settings as _ls
                ag = (_ls() or {}).get("agent", {}) or {}
                tokens = int(ag.get("run_budget_tokens", 0) or 0)
            except Exception:
                tokens = 0
        return max(0, tokens)

    def _run_budget_subagents(self) -> int:
        """Sub-agent ceiling for this run — 0 disables it. OPTIONAL,
        default OFF: absent from the settings and unset on the engine it
        adds no dimension to :meth:`_run_budget_detail` at all.

        **Cost in prompt text:** when it is in force, the wind-down /
        exhausted budget block and the >110 % refusal name it like the
        other dials — one short setting name in text that is already
        injected every turn past 80 % of any configured budget; with the
        cap off the prompt text is unchanged.

        Delegates count towards the PARENT's cap: every sub-agent this
        session ran is in the delegated-spend ledger's session bucket,
        whoever asked for it and whenever it returned. Same precedence as
        the other dimensions: an explicit ``run_budget_subagents``
        attribute beats the settings value ``agent.run_budget_subagents``.
        No command-line flag and no unattended default: off until
        somebody asks for it.
        """
        engine = self._engine
        sub = int(getattr(engine, "run_budget_subagents", 0) or 0)
        if sub <= 0:
            try:
                from delfin.user_settings import load_settings as _ls
                ag = (_ls() or {}).get("agent", {}) or {}
                sub = int(ag.get("run_budget_subagents", 0) or 0)
            except Exception:
                sub = 0
        return max(0, sub)

    def _requests_taken(self) -> int:
        """Model requests this session has actually made.

        Read from the per-session request metrics, one entry per finished
        request (``record_request`` runs once per round, on every
        backend). Requests with a ``purpose`` are side calls (the action
        classifier's verdicts), not the agent's own model requests, and
        are left out — the same split ``format_requests`` makes.
        Best-effort: no metrics file reads as zero, never an exception.
        """
        engine = self._engine
        try:
            from . import turn_metrics as _tm
            session = str(getattr(engine, "session_id", "") or "")
            if not session:
                return 0
            entries = _tm.read_requests(session)
            return sum(1 for e in entries
                       if isinstance(e, dict) and not e.get("purpose"))
        except Exception:
            return 0

    def _run_budget_tokens_spent(self) -> int:
        """Tokens this run has actually spent: the parent's own
        session-cumulative counter PLUS every delegated run's tokens.

        The delegated half comes from the delegated-spend ledger's
        session bucket, the same figure the cost line reports; counting
        only the parent would let a run that delegates never reach its
        token cap. Missing counters (an engine built bare) read as zero.
        """
        engine = self._engine
        usage = getattr(engine, "token_usage", None) or {}
        parent = (int(usage.get("input", 0) or 0)
                  + int(usage.get("output", 0) or 0))
        delegated = 0
        try:
            led = self.delegate_spend().session
            delegated = int(led.input_tokens) + int(led.output_tokens)
        except Exception:
            delegated = 0
        return max(0, parent) + max(0, delegated)

    def _subagents_taken(self) -> int:
        """Sub-agents this session has actually run — the delegated-spend
        ledger's session count, whoever the delegate's parent was.
        """
        try:
            return int(self.delegate_spend().session.count or 0)
        except Exception:
            return 0

    def _run_budget_status(self) -> tuple[float, bool]:
        """(worst fraction spent, exhausted). Fraction is max over the
        enabled dimensions; (0.0, False) when no budget is configured.

        The USD term is a sum over measurable turns only — see
        :meth:`_usd_budget_enforced` for what it cannot see.
        """
        frac, exhausted, _dim = self._run_budget_detail()
        return frac, exhausted

    def _run_budget_detail(self) -> tuple[float, bool, str]:
        """As :meth:`_run_budget_status`, plus the setting that is worst.

        A refusal that names ``run_budget_usd`` when the wall clock or the
        turn count stopped the run sends its reader to the wrong dial.
        """
        engine = self._engine
        usd, secs = self._run_budget()
        turns_cap = self._run_budget_turns()
        frac, dim = 0.0, ""
        if usd > 0:
            share = float(engine.cost_usd) / usd
            if share > frac:
                frac, dim = share, "agent.run_budget_usd"
        if secs > 0:
            import time as _t
            started = float(getattr(engine, "_run_started_at", 0.0) or 0.0)
            if started > 0:
                share = (_t.time() - started) / secs
                if share > frac:
                    frac, dim = share, "agent.run_budget_s"
        if turns_cap > 0:
            share = float(self._turns_taken()) / turns_cap
            if share > frac:
                frac, dim = share, "agent.run_budget_turns"
        reqs_cap = self._run_budget_requests()
        if reqs_cap > 0:
            share = float(self._requests_taken()) / reqs_cap
            if share > frac:
                frac, dim = share, "agent.run_budget_requests"
        tokens_cap = self._run_budget_tokens()
        if tokens_cap > 0:
            share = float(self._run_budget_tokens_spent()) / tokens_cap
            if share > frac:
                frac, dim = share, "agent.run_budget_tokens"
        sub_cap = self._run_budget_subagents()
        if sub_cap > 0:
            share = float(self._subagents_taken()) / sub_cap
            if share > frac:
                frac, dim = share, "agent.run_budget_subagents"
        return frac, frac >= 1.0, dim

    # -- what the model and the operator hear ------------------------------

    def _emit_budget_attention(self, level: str, title: str,
                               detail: str) -> None:
        """Raise ONE durable budget event per level for this run.

        ``budget_warning`` was declared, labelled and ranked in the
        attention inbox with nothing anywhere that emitted it. Everything
        the budget machinery says goes either into the system prompt (the
        model reads it, nobody else) or into the chat stream (gone with
        the scrollback), so an unattended run that stopped on its ceiling
        overnight showed a silent desktop, an empty inbox and a green
        doctor. Best-effort: a notice must never break the turn it
        describes.
        """
        engine = self._engine
        try:
            if level in engine._budget_attention_levels:
                return
            engine._budget_attention_levels.add(level)
            from delfin.agent.attention import emit_attention
            emit_attention(
                "budget_warning",
                session_id=str(getattr(engine, "session_id", "") or ""),
                title=title[:200],
                detail=detail[:400],
                workspace=str(getattr(engine, "repo_dir", "") or ""),
            )
        except Exception:
            pass

    def _build_budget_block(self) -> str:
        """Per-turn budget status for the prompt: silent when unconfigured
        or comfortably below the wind-down threshold; from 80% it instructs
        the agent to wrap up gracefully instead of being cut off mid-work.

        The unmeasured-cost notice comes first and independently of the
        threshold, because the threshold is exactly what an unpriced run
        can never reach: its fraction stays at 0% however long it runs.
        """
        engine = self._engine
        parts: list[str] = []
        try:
            parts.append(engine._unmeasured_budget_block())
        except Exception:
            pass
        try:
            frac, exhausted, dim = self._run_budget_detail()
        except Exception:
            frac, exhausted, dim = 0.0, False, ""
        if frac >= 0.8:
            pct = min(999, int(frac * 100))
            usd, secs = self._run_budget()
            turn_cap = self._run_budget_turns()
            # Every dial in force, so the notice says which run this was
            # and the next line can say which one to turn.
            ceiling = ", ".join(
                [p for p in (f"${usd:.2f}" if usd > 0 else "",
                             f"{secs:.0f}s" if secs > 0 else "",
                             f"{turn_cap} turns" if turn_cap > 0 else "",
                             f"{self._run_budget_requests()} requests"
                             if self._run_budget_requests() > 0 else "",
                             f"{self._run_budget_tokens()} tokens"
                             if self._run_budget_tokens() > 0 else "",
                             f"{self._run_budget_subagents()} sub-agents"
                             if self._run_budget_subagents() > 0 else "")
                 if p])
            # The dial that is actually binding, not a fixed pair: a run
            # stopped by the request cap is not helped by more money.
            raise_it = (f"Raise {dim}" if dim
                        else "Raise agent.run_budget_usd / run_budget_s")
            if exhausted:
                parts.append(
                    "# Run budget EXHAUSTED\n"
                    f"- {pct}% of this run's budget is spent. Finish your "
                    "CURRENT sentence of work only: save state, mark task "
                    "statuses honestly, and summarise what remains. No new "
                    "work."
                )
                engine._emit_budget_attention(
                    "exhausted",
                    "Run budget exhausted — the agent is wrapping up",
                    f"{pct}% of the run budget ({ceiling}) is spent. The "
                    f"agent finishes the current step and starts no new "
                    f"work. {raise_it} or start a new run to continue.")
            else:
                parts.append(
                    "# Run budget wind-down (auto-injected)\n"
                    f"- {pct}% of this run's budget is spent. Start wrapping "
                    "up NOW: complete or checkpoint the current step, "
                    "commit/save state, mark tasks, and summarise open work. "
                    "If a follow-up run makes sense, note exactly where to "
                    "resume. Do not start new large work items."
                )
                engine._emit_budget_attention(
                    "wind_down",
                    "Run budget at 80% — the agent is winding down",
                    f"{pct}% of the run budget ({ceiling}) is spent. The "
                    f"agent is checkpointing and will stop soon.")
        return "\n\n".join(p for p in parts if p)

    def _cost_hard_cap(self) -> float:
        """Per-turn hard cost ceiling in USD -- a runaway circuit-breaker, NOT a
        cumulative session budget. The default is deliberately very high so it
        only ever trips on a loop gone wrong, never on normal use. 0 disables
        it. Configurable via settings ``agent.cost_hard_limit_usd``."""
        try:
            from delfin.user_settings import load_settings as _ls
            v = ((_ls() or {}).get("agent", {}) or {}).get(
                "cost_hard_limit_usd", 50.0)
            return max(0.0, float(v))
        except Exception:
            return 50.0

    # -- what the session delegated ----------------------------------------

    def delegate_spend(self) -> "_agent_metrics.DelegateLedger":
        """This session's delegated spend: per turn, per session, and the
        background part that belongs to no turn.

        Handed out as the ledger rather than as three numbers so a caller
        cannot read the per-turn bucket where it meant the session's.
        """
        engine = self._engine
        ledger = getattr(engine, "_delegate_spend", None)
        if not isinstance(ledger, _agent_metrics.DelegateLedger):
            ledger = _agent_metrics.DelegateLedger()
            engine._delegate_spend = ledger
        return ledger

    def _meter_delegate_spend(self) -> bool:
        """Route every delegated run through a counter that knows which
        turn asked for it.

        The delegate's own telemetry file records tokens and duration but
        carries no parent, session or turn field, so nothing in that file
        can say which turn a delegate belonged to and a timestamp window
        over it would be a guess -- worse, a shared home directory means
        the guess would also pick up another session's delegates. The
        runner the tool executor calls is the one place where a delegate
        and the turn that started it are the same call, so the join is
        made HERE, in memory, and the file is not consulted at all.

        Re-installed at the start of every turn because the client
        re-binds ``subagent_runner`` on every ``set_permissions``. It
        unwraps an existing meter before wrapping, so re-installing can
        never stack two counters over one runner and bill each delegate
        twice.

        Returns True when a meter is in place. False means this backend
        has no permissions object or no runner, i.e. it cannot delegate.
        """
        engine = self._engine
        perms = engine.kit_permissions
        if perms is None:
            return False
        runner = getattr(perms, "subagent_runner", None)
        if runner is None:
            return False
        inner = getattr(runner, "_delegate_meter_inner", None) or runner

        def _metered(*args, **kwargs):
            # Read at ENTRY, not at return: a background delegate comes
            # back long after the turn that started it, and the turn
            # running by then is exactly the one it must NOT be charged
            # to.
            started = engine._turn_id if engine._turn_in_flight else 0
            payload = inner(*args, **kwargs)
            try:
                engine._book_delegate_spend(payload, started_turn=started)
            except Exception:
                pass
            return payload

        _metered._delegate_meter_inner = inner
        try:
            perms.subagent_runner = _metered
        except Exception:
            return False
        return True

    def _book_delegate_spend(self, payload, *, started_turn: int) -> None:
        """Charge one finished delegate, to exactly one place.

        The turn is identified by the pair (a turn is in flight, its id),
        read at the runner's entry and again at its return. Turn ids
        strictly increase and the turn gate allows only one turn in flight
        at a time, so an unchanged pair proves the SAME turn has been
        running throughout -- which is what makes a foreground delegate
        this turn's, and a background one, which returns after its turn
        ended, nobody's.

        Every delegate is booked into the session total and into exactly
        one of the per-turn or the background bucket, so no delegate can
        be counted twice however often the numbers are read afterwards.
        """
        engine = self._engine
        ledger = engine.delegate_spend()
        one = _agent_metrics.DelegateSpend()
        if not one.add_payload(payload,
                               provider=str(getattr(engine, "provider", "")
                                            or "")):
            return
        with engine._lock:
            same_turn = bool(
                started_turn
                and getattr(engine, "_turn_in_flight", False)
                and getattr(engine, "_turn_id", 0) == started_turn)
            (ledger.turn if same_turn else ledger.background).merge(one)
            ledger.session.merge(one)
