from snakebids.snakenull import normalize_inputs_with_snakenull


class DummyRec:
    def __init__(self, entities):
        self.entities = entities


class DummyComponent:
    def __init__(self, requested_wildcards, records):
        self.requested_wildcards = list(requested_wildcards)
        self.matched_files = [DummyRec(e) for e in records]
        self.entities = {}
        self.wildcards = {}


def _entities(component):
    return component.entities


def test_snakenull_disabled_is_noop():
    inputs = {
        "t1w": DummyComponent(
            ["subject", "session", "acquisition"],
            [
                {"subject": "01", "session": "01", "acquisition": "MPRAGE"},
                {"subject": "02"},  # missing session & acquisition
            ],
        )
    }
    cfg = {
        "pybids_inputs": {"t1w": {"wildcards": ["subject", "session", "acquisition"]}},
        "snakenull": {"enabled": False},
    }
    normalize_inputs_with_snakenull(inputs, config=cfg)
    ents = _entities(inputs["t1w"])
    # No processing happened; in particular, no 'snakenull' is injected
    flat = {v for vals in ents.values() for v in vals} if ents else set()
    assert "snakenull" not in flat


def test_snakenull_enables_mixed_entity_normalization():
    inputs = {
        "t1w": DummyComponent(
            ["subject", "session", "acquisition"],
            [
                {"subject": "01", "session": "01", "acquisition": "MPRAGE"},
                {"subject": "02"},  # missing session & acquisition
            ],
        )
    }
    cfg = {
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


def test_snakenull_skips_completely_absent_entities():
    inputs = {
        "t1w": DummyComponent(
            ["subject", "session", "acquisition", "run"],
            [
                {"subject": "01", "session": "01", "acquisition": "MPRAGE"},
                {"subject": "02"},  # still no 'run' anywhere
            ],
        )
    }
    cfg = {
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
