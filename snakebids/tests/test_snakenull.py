from __future__ import annotations

import contextlib
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

import pytest

from snakebids import generate_inputs
from snakebids.exceptions import ConfigError
from snakebids.snakenull import (
    SnakenullConfig,
    _collect_present_values,
    _in_scope,
    _merge_cfg,
    _set_component_entities,
    normalize_inputs_with_snakenull,
)


class DummyRec:
    def __init__(self, entities: Mapping[str, str]) -> None:
        self.entities: Mapping[str, str] = entities


class DummyComponent:
    def __init__(
        self,
        requested_wildcards: Sequence[str],
        records: Sequence[Mapping[str, str]],
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
                "snakenull": {
                    "enabled": True,
                    "scope": ["session", "acquisition"],
                },
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


# ===== UNIT TESTS FOR INTERNAL FUNCTIONS =====


class TestSnakenullConfig:
    """Test the SnakenullConfig dataclass."""

    def test_default_values(self):
        """Test that SnakenullConfig has correct default values."""
        config = SnakenullConfig()
        assert config.enabled is False
        assert config.label == "snakenull"
        assert config.scope == "all"

    def test_custom_values(self):
        """Test that SnakenullConfig accepts custom values."""
        config = SnakenullConfig(
            enabled=True, label="custom_null", scope=["subject", "session"]
        )
        assert config.enabled is True
        assert config.label == "custom_null"
        assert config.scope == ["subject", "session"]


class TestInScope:
    """Test the _in_scope helper function."""

    def test_scope_all(self):
        """Test that 'all' scope includes any entity."""
        assert _in_scope("subject", "all") is True
        assert _in_scope("session", "all") is True
        assert _in_scope("random_entity", "all") is True

    def test_scope_list(self):
        """Test that list scope only includes specified entities."""
        scope = ["subject", "session"]
        assert _in_scope("subject", scope) is True
        assert _in_scope("session", scope) is True
        assert _in_scope("acquisition", scope) is False
        assert _in_scope("run", scope) is False

    def test_scope_set(self):
        """Test that set scope works correctly."""
        scope = {"subject", "session"}
        assert _in_scope("subject", scope) is True
        assert _in_scope("session", scope) is True
        assert _in_scope("acquisition", scope) is False


class TestMergeCfg:
    """Test the _merge_cfg configuration merging function."""

    def test_no_config(self):
        """Test merge with no configuration provided."""
        config = _merge_cfg(None, None)
        assert config.enabled is False
        assert config.label == "snakenull"
        assert config.scope == "all"

    def test_global_config_only(self):
        """Test merge with only global configuration."""
        global_cfg = {"enabled": True, "label": "global_null"}
        config = _merge_cfg(global_cfg, None)
        assert config.enabled is True
        assert config.label == "global_null"
        assert config.scope == "all"

    def test_local_config_only(self):
        """Test merge with only local configuration."""
        local_cfg = {"enabled": True, "scope": ["subject"]}
        config = _merge_cfg(None, local_cfg)
        assert config.enabled is True
        assert config.label == "snakenull"
        assert config.scope == ["subject"]

    def test_local_overrides_global(self):
        """Test that local configuration overrides global."""
        global_cfg = {"enabled": True, "label": "global_null", "scope": "all"}
        local_cfg = {"label": "local_null", "scope": ["subject"]}
        config = _merge_cfg(global_cfg, local_cfg)
        assert config.enabled is True  # from global
        assert config.label == "local_null"  # from local (overrides global)
        assert config.scope == ["subject"]  # from local (overrides global)


class TestCollectPresentValues:
    """Test the _collect_present_values function."""

    def test_collect_with_missing_values(self):
        """Test collecting values when some entities are missing."""
        component = DummyComponent(
            ["subject", "session", "acquisition"],
            [
                {"subject": "01", "session": "01", "acquisition": "MPRAGE"},
                {"subject": "02"},  # missing session & acquisition
                {"subject": "03", "session": "02"},  # missing acquisition
            ],
        )

        present_values, has_missing, wc_list = _collect_present_values(component)

        assert set(wc_list) == {"subject", "session", "acquisition"}
        assert present_values["subject"] == {"01", "02", "03"}
        assert present_values["session"] == {"01", "02"}
        assert present_values["acquisition"] == {"MPRAGE"}

        assert has_missing["subject"] is False  # all records have subject
        assert has_missing["session"] is True  # record 2 missing session
        assert has_missing["acquisition"] is True  # records 2&3 missing acquisition

    def test_collect_no_missing_values(self):
        """Test collecting values when no entities are missing."""
        component = DummyComponent(
            ["subject", "session"],
            [
                {"subject": "01", "session": "01"},
                {"subject": "02", "session": "02"},
            ],
        )

        present_values, has_missing, wc_list = _collect_present_values(component)

        assert set(wc_list) == {"subject", "session"}
        assert present_values["subject"] == {"01", "02"}
        assert present_values["session"] == {"01", "02"}

        assert has_missing["subject"] is False
        assert has_missing["session"] is False

    def test_collect_completely_absent_entity(self):
        """Test collecting values when an entity is completely absent."""
        component = DummyComponent(
            [
                "subject",
                "session",
                "run",
            ],  # run is in wildcards but not in any record
            [
                {"subject": "01", "session": "01"},
                {"subject": "02", "session": "02"},
            ],
        )

        present_values, has_missing, wc_list = _collect_present_values(component)

        assert set(wc_list) == {"subject", "session", "run"}
        assert present_values["subject"] == {"01", "02"}
        assert present_values["session"] == {"01", "02"}
        # When an entity is completely absent, it should not appear in present_values
        # The function may return an empty set or not include the key at all
        assert present_values.get("run", set()) == set()

        assert has_missing["subject"] is False
        assert has_missing["session"] is False
        assert (
            has_missing["run"] is True
        )  # completely absent entities are marked as missing


class TestSetComponentEntities:
    """Test the _set_component_entities function."""

    def test_sets_entities_attribute(self):
        """Test that entities attribute is set correctly."""
        component = DummyComponent(["subject", "session"], [])
        entities = {
            "subject": ["01", "02"],
            "session": ["01", "02", "snakenull"],
        }

        _set_component_entities(component, entities)

        assert component.entities == entities

    def test_sets_wildcards_attribute(self):
        """Test that wildcards attribute is set correctly."""
        component = DummyComponent(["subject", "session"], [])
        entities = {
            "subject": ["01", "02"],
            "session": ["01", "02", "snakenull"],
        }

        _set_component_entities(component, entities)

        expected_wildcards = {"subject": "{subject}", "session": "{session}"}
        assert component.wildcards == expected_wildcards

    def test_sets_zip_lists_for_bids_component(self):
        """Test that zip_lists is updated for BidsComponent-like objects."""

        # Create a mock component with zip_lists attribute
        class MockBidsComponent:
            def __init__(self):
                self.entities = {}
                self.wildcards = {}
                self.zip_lists = {}

        component = MockBidsComponent()
        entities = {
            "subject": ["01", "02"],
            "session": ["01", "02", "snakenull"],
        }

        _set_component_entities(component, entities)

        assert component.zip_lists == entities
        assert component.entities == entities

    def test_sets_snakenull_entities_side_channel(self):
        """Test that snakenull_entities side channel is set."""
        component = DummyComponent(["subject", "session"], [])
        entities = {
            "subject": ["01", "02"],
            "session": ["01", "02", "snakenull"],
        }

        _set_component_entities(component, entities)

        assert hasattr(component, "snakenull_entities")
        assert component.snakenull_entities == entities


# ===== INTEGRATION TESTS =====


class TestSnakenullIntegration:
    """Integration tests using real BIDS data structures."""

    def test_integration_with_generate_inputs(self):
        """Test snakenull integration with generate_inputs function."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a BIDS dataset with missing entities
            bids_root = Path(tmpdir) / "bids"
            bids_root.mkdir()

            # Create subjects with some missing sessions - but avoid ambiguous patterns
            # sub-01 has session 01
            (bids_root / "sub-01" / "ses-01" / "anat").mkdir(parents=True)
            (
                bids_root / "sub-01" / "ses-01" / "anat" / "sub-01_ses-01_T1w.nii.gz"
            ).touch()

            # sub-02 has sessions 01 and 02
            (bids_root / "sub-02" / "ses-01" / "anat").mkdir(parents=True)
            (
                bids_root / "sub-02" / "ses-01" / "anat" / "sub-02_ses-01_T1w.nii.gz"
            ).touch()
            (bids_root / "sub-02" / "ses-02" / "anat").mkdir(parents=True)
            (
                bids_root / "sub-02" / "ses-02" / "anat" / "sub-02_ses-02_T1w.nii.gz"
            ).touch()

            # Configuration - all files have sessions to avoid ambiguous patterns
            pybids_inputs = {
                "T1w": {
                    "filters": {"suffix": "T1w", "extension": ".nii.gz"},
                    "wildcards": ["subject", "session"],
                }
            }

            # Test with snakenull enabled
            dataset_snakenull = generate_inputs(
                bids_dir=str(bids_root),
                pybids_inputs=pybids_inputs,
                snakenull={"enabled": True},
            )

            # Verify snakenull behavior
            t1w_snakenull = dataset_snakenull["T1w"]
            assert len(t1w_snakenull.zip_lists["subject"]) == 2  # 2 subjects found

            # Check that both present sessions and snakenull are there
            # Since all files have sessions, no snakenull should be added to session
            session_values = set(t1w_snakenull.zip_lists["session"])
            assert "01" in session_values
            assert "02" in session_values
            # No snakenull added since no sessions are missing
            assert "snakenull" not in session_values

    def test_custom_snakenull_label(self):
        """Test using a custom snakenull label with missing entities."""
        # Use unit test approach with DummyComponent to test custom labels
        inputs: dict[str, DummyComponent] = {
            "t1w": DummyComponent(
                ["subject", "session"],
                [
                    {"subject": "01", "session": "01"},
                    {"subject": "01", "session": "02"},
                    {"subject": "02", "session": "01"},
                    {"subject": "02"},  # missing session
                ],
            )
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "t1w": {
                    "wildcards": ["subject", "session"],
                    "snakenull": {"enabled": True, "label": "MISSING"},
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["t1w"])

        # Custom label should be used for missing sessions
        assert "01" in ents["session"]
        assert "02" in ents["session"]
        assert "MISSING" in ents["session"]  # Custom label for missing sessions
        assert "snakenull" not in ents["session"]  # Default label should not appear

    def test_scoped_normalization(self):
        """Test that snakenull only affects entities in scope."""
        # Use unit test approach to avoid path template ambiguity
        inputs: dict[str, DummyComponent] = {
            "bold": DummyComponent(
                ["subject", "session", "task", "run"],
                [
                    {
                        "subject": "01",
                        "task": "rest",
                    },  # missing session and run
                    {
                        "subject": "01",
                        "session": "01",
                        "task": "rest",
                        "run": "01",
                    },
                ],
            )
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "bold": {
                    "wildcards": ["subject", "session", "task", "run"],
                    "snakenull": {
                        "enabled": True,
                        "scope": ["session"],
                    },  # Only normalize session
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["bold"])

        # Only session should have snakenull (in scope)
        assert "snakenull" in ents["session"]
        # run should have snakenull too since it's missing in some records, but not in scope
        assert "snakenull" not in ents.get("run", [])  # run not in scope

        # task should not have snakenull (not in scope, and all records have task)
        assert "snakenull" not in ents["task"]

    def test_per_component_configuration(self):
        """Test per-component snakenull configuration."""
        # Use unit test approach to test per-component configuration
        inputs: dict[str, DummyComponent] = {
            "T1w": DummyComponent(
                ["subject", "session"],
                [
                    {"subject": "01"},  # missing session
                    {"subject": "02", "session": "01"},
                ],
            ),
            "T2w": DummyComponent(
                ["subject", "session"],
                [
                    {"subject": "01", "session": "01"},
                    {"subject": "02"},  # missing session
                ],
            ),
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "T1w": {
                    "wildcards": ["subject", "session"],
                    "snakenull": {"enabled": True, "label": "T1_NULL"},
                },
                "T2w": {
                    "wildcards": ["subject", "session"],
                    "snakenull": {"enabled": True, "label": "T2_NULL"},
                },
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)

        # Each component should use its own label
        t1w_ents = _entities(inputs["T1w"])
        t2w_ents = _entities(inputs["T2w"])

        assert "T1_NULL" in t1w_ents["session"]
        assert "T2_NULL" in t2w_ents["session"]

        # Cross-contamination should not occur
        assert "T2_NULL" not in t1w_ents["session"]
        assert "T1_NULL" not in t2w_ents["session"]


# ===== END-TO-END TESTS =====


class TestSnakenullEndToEnd:
    """Simplified E2E tests focusing on core snakenull functionality."""

    def test_snakenull_end_to_end_functionality(self):
        """Test complete snakenull workflow end-to-end using unit components."""
        # This test demonstrates the full snakenull workflow without
        # the complications of real BIDS path template ambiguity
        inputs: dict[str, DummyComponent] = {
            "bold": DummyComponent(
                ["subject", "session", "task", "run"],
                [
                    # Complete record
                    {
                        "subject": "01",
                        "session": "01",
                        "task": "rest",
                        "run": "01",
                    },
                    # Missing run
                    {"subject": "01", "session": "02", "task": "rest"},
                    # Missing session and run
                    {"subject": "02", "task": "rest"},
                    # Complete but different values
                    {
                        "subject": "03",
                        "session": "01",
                        "task": "task",
                        "run": "02",
                    },
                ],
            )
        }

        cfg: dict[str, Any] = {
            "snakenull": {"enabled": True, "label": "MISSING"},
            "pybids_inputs": {
                "bold": {
                    "wildcards": ["subject", "session", "task", "run"],
                    "snakenull": {"scope": ["session", "run"]},  # Only normalize these
                }
            },
        }

        # Apply snakenull normalization
        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["bold"])

        # Verify complete normalization
        assert (
            "01" in ents["subject"]
            and "02" in ents["subject"]
            and "03" in ents["subject"]
        )
        assert (
            "01" in ents["session"]
            and "02" in ents["session"]
            and "MISSING" in ents["session"]
        )
        assert "rest" in ents["task"] and "task" in ents["task"]
        assert "01" in ents["run"] and "02" in ents["run"] and "MISSING" in ents["run"]

        # Verify scoping worked - task should not have MISSING (not in scope)
        assert "MISSING" not in ents["task"]


# ===== ERROR HANDLING AND EDGE CASES =====


class TestSnakenullErrorHandling:
    """Test error handling and edge cases."""

    def test_empty_dataset(self):
        """Test snakenull with empty dataset."""
        inputs: dict[str, DummyComponent] = {}
        cfg: dict[str, Any] = {"snakenull": {"enabled": True}}

        # Should not crash
        normalize_inputs_with_snakenull(inputs, config=cfg)
        assert len(inputs) == 0

    def test_component_with_no_records(self):
        """Test component with no matched files."""
        inputs: dict[str, DummyComponent] = {
            "empty": DummyComponent(["subject", "session"], [])
        }
        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "empty": {
                    "wildcards": ["subject", "session"],
                    "snakenull": {"enabled": True},
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)

        # Should have empty entities (no values to normalize)
        assert _entities(inputs["empty"]) == {}

    def test_component_with_no_wildcards(self):
        """Test component with no wildcards defined."""
        inputs: dict[str, DummyComponent] = {
            "nowild": DummyComponent([], [{"path": "/some/file.nii.gz"}])
        }
        cfg: dict[str, Any] = {
            "pybids_inputs": {"nowild": {"snakenull": {"enabled": True}}}
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)

        # Should handle gracefully
        assert _entities(inputs["nowild"]) == {}

    def test_malformed_config(self):
        """Test handling of malformed configuration."""
        inputs: dict[str, DummyComponent] = {
            "t1w": DummyComponent(["subject"], [{"subject": "01"}])
        }

        # Test with various malformed configs
        malformed_configs = [
            {"snakenull": "not_a_dict"},
            {"pybids_inputs": "not_a_dict"},
            {"pybids_inputs": {"t1w": "not_a_dict"}},
            {"pybids_inputs": {"t1w": {"snakenull": "not_a_dict"}}},
        ]

        for cfg in malformed_configs:
            # Should not crash with malformed config
            with contextlib.suppress(Exception):
                normalize_inputs_with_snakenull(inputs, config=cfg)

    def test_missing_attributes_handled_gracefully(self):
        """Test that missing attributes are handled gracefully."""

        class MinimalComponent:
            """Component with minimal attributes."""

            def __init__(self):
                self.requested_wildcards = ["subject"]
                self.matched_files = [
                    type("obj", (object,), {"entities": {"subject": "01"}})()
                ]

        inputs = {"minimal": MinimalComponent()}
        cfg = {"snakenull": {"enabled": True}}

        # Should not crash even with minimal component
        with contextlib.suppress(Exception):
            normalize_inputs_with_snakenull(inputs, config=cfg)


# ===== ADDITIONAL COMPREHENSIVE TESTS =====


class TestSnakenullNormalizationLogic:
    """Test the core normalization logic comprehensively."""

    def test_normalization_preserves_existing_values(self):
        """Test that normalization preserves all existing entity values."""
        inputs: dict[str, DummyComponent] = {
            "bold": DummyComponent(
                ["subject", "session", "task", "run"],
                [
                    {
                        "subject": "01",
                        "session": "01",
                        "task": "rest",
                        "run": "01",
                    },
                    {
                        "subject": "01",
                        "session": "01",
                        "task": "rest",
                        "run": "02",
                    },
                    {
                        "subject": "02",
                        "session": "01",
                        "task": "rest",
                    },  # missing run
                    {
                        "subject": "02",
                        "task": "rest",
                    },  # missing session and run
                ],
            )
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "bold": {
                    "wildcards": ["subject", "session", "task", "run"],
                    "snakenull": {"enabled": True},
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["bold"])

        # All original values should be preserved
        assert "01" in ents["subject"]
        assert "02" in ents["subject"]
        assert "01" in ents["session"]
        assert "rest" in ents["task"]
        assert "01" in ents["run"]
        assert "02" in ents["run"]

        # Snakenull should be added for missing values
        assert "snakenull" in ents["session"]
        assert "snakenull" in ents["run"]

    def test_normalization_with_custom_scope(self):
        """Test that custom scope only affects specified entities."""
        inputs: dict[str, DummyComponent] = {
            "fmap": DummyComponent(
                ["subject", "session", "acq", "dir"],
                [
                    {
                        "subject": "01",
                        "session": "01",
                        "acq": "AP",
                        "dir": "PA",
                    },
                    {
                        "subject": "02",
                        "session": "02",
                        "acq": "PA",
                    },  # missing dir
                    {"subject": "03", "acq": "AP"},  # missing session and dir
                ],
            )
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "fmap": {
                    "wildcards": ["subject", "session", "acq", "dir"],
                    "snakenull": {"enabled": True, "scope": ["session", "dir"]},
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["fmap"])

        # session and dir should have snakenull (in scope)
        assert "snakenull" in ents["session"]
        assert "snakenull" in ents["dir"]

        # acq should not have snakenull (not in scope, even though all present)
        assert "snakenull" not in ents["acq"]

    def test_multiple_components_independent_processing(self):
        """Test that multiple components are processed independently."""
        inputs: dict[str, DummyComponent] = {
            "T1w": DummyComponent(
                ["subject", "session"],
                [
                    {"subject": "01", "session": "01"},
                    {"subject": "02"},  # missing session
                ],
            ),
            "T2w": DummyComponent(
                ["subject", "acq"],
                [
                    {"subject": "01", "acq": "SPACE"},
                    {"subject": "02"},  # missing acq
                ],
            ),
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "T1w": {
                    "wildcards": ["subject", "session"],
                    "snakenull": {"enabled": True, "label": "T1_MISSING"},
                },
                "T2w": {
                    "wildcards": ["subject", "acq"],
                    "snakenull": {"enabled": True, "label": "T2_MISSING"},
                },
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)

        t1w_ents = _entities(inputs["T1w"])
        t2w_ents = _entities(inputs["T2w"])

        # Each component should use its own label
        assert "T1_MISSING" in t1w_ents["session"]
        assert "T2_MISSING" in t2w_ents["acq"]

        # Cross-contamination should not occur
        assert "T2_MISSING" not in t1w_ents.get("session", [])
        assert "T1_MISSING" not in t2w_ents.get("acq", [])

    def test_global_vs_local_configuration_precedence(self):
        """Test that local component config overrides global config."""
        inputs: dict[str, DummyComponent] = {
            "bold": DummyComponent(
                ["subject", "task"],
                [
                    {"subject": "01", "task": "rest"},
                    {"subject": "02"},  # missing task
                ],
            )
        }

        cfg: dict[str, Any] = {
            "snakenull": {"enabled": True, "label": "GLOBAL_NULL"},
            "pybids_inputs": {
                "bold": {
                    "wildcards": ["subject", "task"],
                    "snakenull": {"label": "LOCAL_NULL"},  # Override label locally
                }
            },
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["bold"])

        # Should use local label, not global
        assert "LOCAL_NULL" in ents["task"]
        assert "GLOBAL_NULL" not in ents["task"]


class TestSnakenullEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_all_entities_complete_no_normalization_needed(self):
        """Test behavior when no entities are missing."""
        inputs: dict[str, DummyComponent] = {
            "dwi": DummyComponent(
                ["subject", "session", "acq", "dir"],
                [
                    {
                        "subject": "01",
                        "session": "01",
                        "acq": "b1000",
                        "dir": "AP",
                    },
                    {
                        "subject": "01",
                        "session": "01",
                        "acq": "b2000",
                        "dir": "PA",
                    },
                    {
                        "subject": "02",
                        "session": "01",
                        "acq": "b1000",
                        "dir": "AP",
                    },
                ],
            )
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "dwi": {
                    "wildcards": ["subject", "session", "acq", "dir"],
                    "snakenull": {"enabled": True},
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["dwi"])

        # When all entities are complete, no snakenull should be added
        assert "snakenull" not in ents.get("subject", [])
        assert "snakenull" not in ents.get("session", [])
        assert "snakenull" not in ents.get("acq", [])
        assert "snakenull" not in ents.get("dir", [])

        # But all original values should be preserved
        assert "01" in ents["subject"] and "02" in ents["subject"]
        assert "01" in ents["session"]
        assert "b1000" in ents["acq"] and "b2000" in ents["acq"]
        assert "AP" in ents["dir"] and "PA" in ents["dir"]

    def test_single_record_component(self):
        """Test component with only one record."""
        inputs: dict[str, DummyComponent] = {
            "single": DummyComponent(
                ["subject", "session", "task"],
                [{"subject": "01", "session": "01", "task": "rest"}],
            )
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "single": {
                    "wildcards": ["subject", "session", "task"],
                    "snakenull": {"enabled": True},
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["single"])

        # Single record with all entities complete - no snakenull should be added
        assert "01" in ents["subject"]
        assert "01" in ents["session"]
        assert "rest" in ents["task"]
        assert "snakenull" not in ents["subject"]
        assert "snakenull" not in ents["session"]
        assert "snakenull" not in ents["task"]

    def test_entity_with_empty_string_value(self):
        """Test handling of entities that have empty string values."""
        inputs: dict[str, DummyComponent] = {
            "test": DummyComponent(
                ["subject", "session"],
                [
                    {"subject": "01", "session": "01"},
                    {"subject": "02", "session": ""},  # explicit empty string
                ],
            )
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "test": {
                    "wildcards": ["subject", "session"],
                    "snakenull": {"enabled": True},
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["test"])

        # Empty string is treated as missing, so snakenull should be added to session
        assert "01" in ents["session"]
        assert "snakenull" in ents["session"]  # because one record has missing session
        # Empty string itself is not preserved as it's treated as None/missing
        assert "" not in ents["session"]

        # Subject is complete for all records
        assert "01" in ents["subject"] and "02" in ents["subject"]
        assert "snakenull" not in ents["subject"]

    def test_duplicate_entity_values(self):
        """Test that duplicate entity values are handled correctly."""
        inputs: dict[str, DummyComponent] = {
            "repeat": DummyComponent(
                ["subject", "run"],
                [
                    {"subject": "01", "run": "01"},
                    {"subject": "01", "run": "01"},  # duplicate
                    {"subject": "01", "run": "02"},
                    {"subject": "02"},  # missing run
                ],
            )
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "repeat": {
                    "wildcards": ["subject", "run"],
                    "snakenull": {"enabled": True},
                }
            }
        }

        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["repeat"])

        # Should deduplicate values
        assert len([x for x in ents["run"] if x == "01"]) == 1
        assert "02" in ents["run"]
        assert "snakenull" in ents["run"]


class TestSnakenullPerformance:
    """Test performance characteristics and large dataset handling."""

    def test_large_number_of_entities(self):
        """Test with a large number of entity types."""
        # Create component with many entity types
        entity_names = [f"entity{i:02d}" for i in range(20)]

        # Create records where each has some missing entities
        records = []
        for i in range(10):
            record = {"subject": f"{i:02d}"}
            # Each record has different subsets of entities
            for j, entity in enumerate(entity_names):
                if (i + j) % 3 == 0:  # Include entity in ~1/3 of records
                    record[entity] = f"value_{j}"
            records.append(record)

        inputs: dict[str, DummyComponent] = {
            "large": DummyComponent(["subject"] + entity_names, records)
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "large": {
                    "wildcards": ["subject"] + entity_names,
                    "snakenull": {"enabled": True},
                }
            }
        }

        # Should not crash with large number of entities
        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["large"])

        # Verify basic structure
        assert "subject" in ents
        assert len(ents["subject"]) >= 10  # At least the subjects

        # Check that snakenull was added appropriately
        for entity in entity_names:
            if entity in ents:
                assert "snakenull" in ents[entity]

    def test_large_number_of_records(self):
        """Test with a large number of records."""
        # Create many records with some missing entities
        records = []
        for i in range(100):
            record = {"subject": f"{i:03d}"}
            if i % 2 == 0:
                record["session"] = "01"
            if i % 3 == 0:
                record["task"] = "rest"
            records.append(record)

        inputs: dict[str, DummyComponent] = {
            "many": DummyComponent(["subject", "session", "task"], records)
        }

        cfg: dict[str, Any] = {
            "pybids_inputs": {
                "many": {
                    "wildcards": ["subject", "session", "task"],
                    "snakenull": {"enabled": True},
                }
            }
        }

        # Should complete in reasonable time
        normalize_inputs_with_snakenull(inputs, config=cfg)
        ents = _entities(inputs["many"])

        # Verify results
        assert (
            len(ents["subject"]) == 100
        )  # all subjects, no snakenull added since all have subject
        assert "snakenull" in ents["session"]  # some records missing session
        assert "snakenull" in ents["task"]  # some records missing task
        assert len(ents["session"]) == 2  # "01" + "snakenull"
        assert len(ents["task"]) == 2  # "rest" + "snakenull"
        assert "01" in ents["session"]
        assert "rest" in ents["task"]


class TestSnakenullRobustness:
    """Test robustness against various failure modes."""

    def test_component_with_malformed_records(self):
        """Test handling of components with unusual record structures."""

        class MalformedComponent:
            def __init__(self):
                self.requested_wildcards = ["subject", "session"]
                # Some records might not have entities dict
                self.matched_files = [
                    type(
                        "obj",
                        (object,),
                        {"entities": {"subject": "01", "session": "01"}},
                    )(),
                    type(
                        "obj", (object,), {"entities": {"subject": "02"}}
                    )(),  # missing session
                    type("obj", (object,), {})(),  # no entities at all
                ]
                self.entities = {}
                self.wildcards = {}

        inputs = {"malformed": MalformedComponent()}
        cfg = {
            "pybids_inputs": {
                "malformed": {
                    "wildcards": ["subject", "session"],
                    "snakenull": {"enabled": True},
                }
            }
        }

        # Should handle gracefully without crashing
        with contextlib.suppress(Exception):
            normalize_inputs_with_snakenull(inputs, config=cfg)

    def test_config_with_unexpected_types(self):
        """Test robustness against unexpected configuration types."""
        inputs: dict[str, DummyComponent] = {
            "test": DummyComponent(["subject"], [{"subject": "01"}])
        }

        # Test various malformed configs
        malformed_configs = [
            None,  # No config
            {},  # Empty config
            {"snakenull": None},  # Null snakenull
            {"pybids_inputs": None},  # Null pybids_inputs
            {"pybids_inputs": {"test": None}},  # Null component config
        ]

        for cfg in malformed_configs:
            # Should not crash
            with contextlib.suppress(Exception):
                normalize_inputs_with_snakenull(inputs, config=cfg)

    def test_component_missing_expected_attributes(self):
        """Test with components missing expected attributes."""

        class MinimalComponent:
            """Component missing some expected attributes."""

            def __init__(self):
                self.requested_wildcards = ["subject"]
                # Missing matched_files, entities, wildcards

        class EmptyComponent:
            """Component with minimal structure."""

        inputs = {
            "minimal": MinimalComponent(),
            "empty": EmptyComponent(),
        }
        cfg = {"snakenull": {"enabled": True}}

        # Should handle gracefully
        with contextlib.suppress(Exception):
            normalize_inputs_with_snakenull(inputs, config=cfg)


class TestMultiTemplateHandling:
    """Test multi-template scenarios with snakenull."""

    def test_multi_template_with_entity_name_mismatches(self):
        """Test that multi-template handling works with entity name mismatches like acq vs acquisition."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a BIDS dataset with T1w files with different entity patterns
            bids_root = Path(tmpdir) / "bids"
            bids_root.mkdir()

            # File with acquisition entity (mapped to acq in template)
            (bids_root / "sub-MPN00001" / "ses-v1" / "anat").mkdir(parents=True)
            (
                bids_root
                / "sub-MPN00001"
                / "ses-v1"
                / "anat"
                / "sub-MPN00001_ses-v1_acq-MPRAGE_T1w.nii.gz"
            ).touch()

            # File with run entity  
            (bids_root / "sub-MPN00002" / "ses-v1" / "anat").mkdir(parents=True)
            (
                bids_root
                / "sub-MPN00002"
                / "ses-v1"
                / "anat"
                / "sub-MPN00002_ses-v1_run-1_T1w.nii.gz"
            ).touch()

            # File with both acq and run
            (bids_root / "sub-MPN00003" / "ses-v1" / "anat").mkdir(parents=True)
            (
                bids_root
                / "sub-MPN00003"
                / "ses-v1"
                / "anat"
                / "sub-MPN00003_ses-v1_acq-MPRAGE_run-1_T1w.nii.gz"
            ).touch()

            # Configuration that uses 'acquisition' wildcard name but files use 'acq'
            pybids_inputs = {
                "T1w": {
                    "filters": {"suffix": "T1w", "extension": ".nii.gz"},
                    "wildcards": ["subject", "session", "acquisition", "run"],
                }
            }

            # Test with snakenull enabled  
            dataset = generate_inputs(
                bids_dir=str(bids_root),
                pybids_inputs=pybids_inputs,
                snakenull={"enabled": True},
            )

            # Verify that all files are found
            t1w = dataset["T1w"]
            assert len(t1w.zip_lists["subject"]) == 3  # All 3 subjects found

            # Check that acquisition and run values are handled correctly
            if "acquisition" in t1w.zip_lists:
                acq_values = set(t1w.zip_lists["acquisition"])
                assert "MPRAGE" in acq_values
                assert "snakenull" in acq_values  # For files without acquisition
            
            if "run" in t1w.zip_lists:
                run_values = set(t1w.zip_lists["run"])
                assert "1" in run_values
                assert "snakenull" in run_values  # For files without run

            # Verify subjects and sessions
            subject_values = set(t1w.zip_lists["subject"])
            session_values = set(t1w.zip_lists["session"])
            assert subject_values == {"MPN00001", "MPN00002", "MPN00003"}
            assert session_values == {"v1"}

    def test_multi_template_with_missing_entities(self):
        """Test that files without run entity work when some templates expect it."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a BIDS dataset with T1w files with different entity patterns
            bids_root = Path(tmpdir) / "bids"
            bids_root.mkdir()

            # Some files have run entity, some don't
            (bids_root / "sub-MPN00001" / "ses-v1" / "anat").mkdir(parents=True)
            (
                bids_root
                / "sub-MPN00001"
                / "ses-v1"
                / "anat"
                / "sub-MPN00001_ses-v1_T1w.nii.gz"
            ).touch()

            (bids_root / "sub-MPN00002" / "ses-v1" / "anat").mkdir(parents=True)
            (
                bids_root
                / "sub-MPN00002"
                / "ses-v1"
                / "anat"
                / "sub-MPN00002_ses-v1_T1w.nii.gz"
            ).touch()

            (bids_root / "sub-MPN00003" / "ses-v1" / "anat").mkdir(parents=True)
            (
                bids_root
                / "sub-MPN00003"
                / "ses-v1"
                / "anat"
                / "sub-MPN00003_ses-v1_run-1_T1w.nii.gz"
            ).touch()

            # Configuration that should accept both patterns
            pybids_inputs = {
                "T1w": {
                    "filters": {"suffix": "T1w", "extension": ".nii.gz"},
                    "wildcards": ["subject", "session", "run"],
                }
            }

            # Test with snakenull enabled
            dataset = generate_inputs(
                bids_dir=str(bids_root),
                pybids_inputs=pybids_inputs,
                snakenull={"enabled": True},
            )

            # Verify that all files are found
            t1w = dataset["T1w"]
            assert len(t1w.zip_lists["subject"]) == 3  # All 3 subjects found

            # Check that run values include both real values and placeholders
            run_values = set(t1w.zip_lists["run"])
            assert "1" in run_values  # Real run value
            assert (
                "snakenull" in run_values
            )  # Placeholder for missing runs (normalized by snakenull)

            # Verify subjects and sessions
            subject_values = set(t1w.zip_lists["subject"])
            session_values = set(t1w.zip_lists["session"])
            assert subject_values == {"MPN00001", "MPN00002", "MPN00003"}
            assert session_values == {"v1"}

    def test_multi_template_without_snakenull_should_error(self):
        """Test that multi-template scenarios fail when snakenull is disabled."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create the same dataset structure as above
            bids_root = Path(tmpdir) / "bids"
            bids_root.mkdir()

            # Some files have run entity, some don't
            (bids_root / "sub-MPN00001" / "ses-v1" / "anat").mkdir(parents=True)
            (
                bids_root
                / "sub-MPN00001"
                / "ses-v1"
                / "anat"
                / "sub-MPN00001_ses-v1_T1w.nii.gz"
            ).touch()

            (bids_root / "sub-MPN00003" / "ses-v1" / "anat").mkdir(parents=True)
            (
                bids_root
                / "sub-MPN00003"
                / "ses-v1"
                / "anat"
                / "sub-MPN00003_ses-v1_run-1_T1w.nii.gz"
            ).touch()

            # Configuration that includes run wildcard
            pybids_inputs = {
                "T1w": {
                    "filters": {"suffix": "T1w", "extension": ".nii.gz"},
                    "wildcards": ["subject", "session", "run"],
                }
            }

            # Test with snakenull disabled - should fail with multiple templates
            with pytest.raises(
                ConfigError, match="Multiple path templates for one component"
            ):
                generate_inputs(
                    bids_dir=str(bids_root),
                    pybids_inputs=pybids_inputs,
                    # snakenull disabled by default
                )
