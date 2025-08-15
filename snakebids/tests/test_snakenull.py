from __future__ import annotations

from typing import Any, Mapping, Sequence

from snakebids.snakenull import normalize_inputs_with_snakenull


class DummyRec:
    def __init__(self, entities: Mapping[str, str]) -> None:
        self.entities: Mapping[str, str] = entities


class DummyComponent:
    def __init__(
        self, requested_wildcards: Sequence[str], records: Sequence[Mapping[str, str]]
    ) -> None:
        self.requested_wildcards: list[str] = list(requested_wildcards)
        self.matched_files: list[DummyRec] = [DummyRec(e) for e in records]
        # Filled by normalize_inputs_with_snakenull; mapping entity -> list of values
        self.entities: dict[str, list[str]] = {}
        # Mapping entity -> "{entity}" placeholder
        self.wildcards: dict[str, str] = {}


def _entities(component: DummyComponent) -> dict[str, list[str]]:
    return component.entities


def test_snakenull_disabled_is_noop() -> None:
    inputs: dict[str, DummyComponent] = {
        "t1w": DummyComponent(
            ["subject", "session", "acquisition"],
            [
                {"subject": "01", "session": "01", "acquisition": "MPRAGE"},
                {"subject": "02"},  # missing session & acquisition
            ],
        )
    }
    cfg: dict[str, Any] = {
        "pybids_inputs": {"t1w": {"wildcards": ["subject", "session", "acquisition"]}},
        "snakenull": {"enabled": False},
    }
    normalize_inputs_with_snakenull(inputs, config=cfg)
    ents = _entities(inputs["t1w"])
    # No processing happened; in particular, no 'snakenull' is injected
    flat: set[str] = {v for vals in ents.values() for v in vals} if ents else set()
    assert "snakenull" not in flat


def test_snakenull_enables_mixed_entity_normalization() -> None:
    inputs: dict[str, DummyComponent] = {
        "t1w": DummyComponent(
            ["subject", "session", "acquisition"],
            [
                {"subject": "01", "session": "01", "acquisition": "MPRAGE"},
                {"subject": "02"},  # missing session & acquisition
            ],
        )
    }
    cfg: dict[str, Any] = {
        "pybids_inputs": {
            "t1w": {
                "wildcards": ["subject", "session", "acquisition"],
                "snakenull": {"enabled": True, "scope": ["session", "acquisition"]},
            }
        }
    }
    normalize_inputs_with_snakenull(inputs, config=cfg)
    ents = _entities(inputs["t1w"])
    assert "acquisition" in ents
    assert set(ents["acquisition"]) == {"MPRAGE", "snakenull"}
    assert "session" in ents
    assert set(ents["session"]) == {"01", "snakenull"}


def test_snakenull_skips_completely_absent_entities() -> None:
    inputs: dict[str, DummyComponent] = {
        "t1w": DummyComponent(
            ["subject", "session", "acquisition", "run"],
            [
                {"subject": "01", "session": "01", "acquisition": "MPRAGE"},
                {"subject": "02"},  # still no 'run' anywhere
            ],
        )
    }
    cfg: dict[str, Any] = {
        "pybids_inputs": {
            "t1w": {
                "wildcards": ["subject", "session", "acquisition", "run"],
                "snakenull": {"enabled": True, "scope": "all"},
            }
        }
    }
    normalize_inputs_with_snakenull(inputs, config=cfg)
    ents = _entities(inputs["t1w"])
    assert "run" not in ents
