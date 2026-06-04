"""Tests for delfin.fffree.grip_loss_weights_tuned (Mission 2 helper)."""
from __future__ import annotations

import os

import pytest

import delfin.fffree.grip_loss_weights_tuned as gwt


def test_env_flag_off_returns_empty(monkeypatch):
    monkeypatch.delenv(gwt.ENV_FLAG, raising=False)
    assert gwt.get_loss_weights() == {}


def test_env_flag_on_returns_table(monkeypatch):
    monkeypatch.setenv(gwt.ENV_FLAG, "1")
    out = gwt.get_loss_weights()
    assert "bond:M-donor" in out
    assert out["bond:M-donor"] == 7.0
    assert out["angle:donor-M-donor"] == 3.0


def test_apply_weights_noop_when_off(monkeypatch):
    monkeypatch.delenv(gwt.ENV_FLAG, raising=False)
    d = {"bond": 4.5}
    out = gwt.apply_weights(d)
    assert out is d
    assert d == {"bond": 4.5}


def test_apply_weights_overrides_class_keys(monkeypatch):
    monkeypatch.setenv(gwt.ENV_FLAG, "1")
    d = {"bond": 4.5, "bond:M-donor": 5.0}
    out = gwt.apply_weights(d)
    assert out["bond:M-donor"] == 7.0
    # default global "bond" stays as it was (caller-provided)
    assert out["bond"] == 4.5


def test_apply_weights_fills_defaults_when_allowed(monkeypatch):
    monkeypatch.setenv(gwt.ENV_FLAG, "1")
    d = {}  # caller provides nothing
    out = gwt.apply_weights(d, allow_default=True)
    # global defaults inserted
    assert out["bond"] == 5.0
    assert out["angle"] == 2.0


def test_apply_weights_no_default_fill_when_disabled(monkeypatch):
    monkeypatch.setenv(gwt.ENV_FLAG, "1")
    d = {}
    out = gwt.apply_weights(d, allow_default=False)
    # only class-specific overrides
    assert "bond" not in out
    assert "bond:M-donor" in out


def test_env_flag_constant():
    assert gwt.ENV_FLAG == "DELFIN_GRIP_LOSS_WEIGHTS_TUNED"


def test_weights_tuned_immutable():
    # mutating a returned copy must not affect the canonical table
    out = gwt.WEIGHTS_TUNED.copy()
    out["bond:M-donor"] = 99.0
    assert gwt.WEIGHTS_TUNED["bond:M-donor"] == 7.0
