"""
Post-processing utilities to normalize mixed/absent entities in Snakebids inputs.

This module assigns a configurable placeholder label (default: "snakenull") to missing
entities and removes entities that are entirely absent in a component. It operates on
the public objects returned by `snakebids.generate_inputs()` and does not require
internal hooks.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping, MutableMapping


@dataclass(frozen=True)
class SnakenullConfig:
    """Configuration for snakenull normalization.

    Attributes
    ----------
    enabled:
        Turn normalization on for this component or globally.
    label:
        Placeholder value assigned to missing entities.
    include_prefix:
        If True, include the key when rendering paths (e.g., ``_acq-snakenull_``);
        if False, omit the entire segment.
    scope:
        Either ``"all"`` or an iterable of entity names to normalize.
    """

    enabled: bool = False
    label: str = "snakenull"
    include_prefix: bool = True  # downstream path builders may use this
    scope: Iterable[str] | str = "all"  # "all" or iterable of entity names


def _in_scope(ent: str, scope: Iterable[str] | str) -> bool:
    """Return True if an entity name is within the configured normalization scope."""
    return scope == "all" or ent in set(scope)


def _merge_cfg(
    global_cfg: Mapping[str, Any] | None, local_cfg: Mapping[str, Any] | None
) -> SnakenullConfig:
    """Merge global and per-component snakenull configuration into a dataclass."""
    base = {
        "enabled": False,
        "label": "snakenull",
        "include_prefix": True,
        "scope": "all",
    }
    if global_cfg:
        base.update(global_cfg)
    if local_cfg:
        base.update(local_cfg)
    return SnakenullConfig(**base)


def _wildcards_list_from_component(component) -> list[str]:
    """Return wildcard names from the component in a tolerant way."""
    if hasattr(component, "requested_wildcards"):
        return list(component.requested_wildcards)
    if hasattr(component, "wildcards") and isinstance(component.wildcards, Mapping):
        return list(component.wildcards.keys())
    if isinstance(component, Mapping) and "wildcards" in component:
        w = component["wildcards"]
        if isinstance(w, Mapping):
            return list(w.keys())
        if isinstance(w, (list, tuple, set)):
            return list(w)
    return []


def _files_from_component(component) -> list:
    """Return the iterable of matched records/files from the component."""
    for attr in ("matched_files", "files", "items", "records"):
        if hasattr(component, attr):
            return list(getattr(component, attr) or [])
    if isinstance(component, Mapping):
        for key in ("matched_files", "files", "items", "records"):
            if key in component:
                return list(component[key] or [])
    return []


def _entities_from_record(rec) -> Mapping[str, Any]:
    """Return the entities mapping from a PyBIDS-like record."""
    if hasattr(rec, "entities"):
        ents = rec.entities
        return ents if isinstance(ents, Mapping) else {}
    if isinstance(rec, Mapping) and "entities" in rec:
        ents = rec.get("entities")
        return ents if isinstance(ents, Mapping) else {}
    return {}


def _collect_present_values(
    component,
) -> tuple[dict[str, set[str]], dict[str, bool], list[str]]:
    """Collect present entity values and detect missing ones across matched files."""
    wc_list = _wildcards_list_from_component(component)
    present_values: dict[str, set[str]] = {w: set() for w in wc_list}
    has_missing: dict[str, bool] = {w: False for w in wc_list}

    for rec in _files_from_component(component):
        ents = _entities_from_record(rec)
        for ent in wc_list:
            val = ents.get(ent)
            if val not in (None, ""):
                present_values[ent].add(str(val))
            else:
                has_missing[ent] = True

    return present_values, has_missing, wc_list


def _set_component_entities(component, entities: Mapping[str, list[str]]) -> None:
    """Write back normalized entities and rebuild wildcards if needed."""
    # entities
    if hasattr(component, "entities"):
        component.entities = dict(entities)
    elif isinstance(component, MutableMapping):
        component["entities"] = dict(entities)
    # wildcards
    wildcards = {k: "{" + k + "}" for k in entities}
    if hasattr(component, "wildcards"):
        component.wildcards = wildcards
    elif isinstance(component, MutableMapping):
        component["wildcards"] = wildcards
    # optional flags for downstream path builders
    if hasattr(component, "__dict__"):
        component.snakenull_label = None
        component.snakenull_include_prefix = True
    if isinstance(component, MutableMapping):
        component["snakenull_label"] = None
        component["snakenull_include_prefix"] = True


def normalize_inputs_with_snakenull(
    inputs: Mapping[str, Any],
    *,
    config: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    """Normalize mixed/absent entities on a Snakebids inputs mapping in-place.

    Post-processes the result of `snakebids.generate_inputs()` to:
      1) Remove entities listed in wildcards that are entirely absent in the dataset.
      2) For entities present in some files but missing in others, insert a placeholder
         value (default: ``"snakenull"``).

    Parameters
    ----------
    inputs
        Object returned by `snakebids.generate_inputs()`, mapping component name to
        component.
    config
        Top-level Snakebids config dict. Read from:
        - `config["snakenull"]` (global defaults)
        - `config["pybids_inputs"][<component>]["snakenull"]` (per-component override)

    Returns
    -------
    Mapping[str, Any]
        The same mapping (mutated in place).
    """
    pybids_inputs = {}
    if config and isinstance(config.get("pybids_inputs"), Mapping):
        pybids_inputs = config["pybids_inputs"]
    global_cfg = {}
    if config and isinstance(config.get("snakenull"), Mapping):
        global_cfg = config["snakenull"]

    for cname, comp in inputs.items():
        local_cfg = {}
        if isinstance(pybids_inputs.get(cname), Mapping):
            local_cfg = pybids_inputs[cname].get("snakenull", {}) or {}
        s_cfg = _merge_cfg(global_cfg, local_cfg)

        if not s_cfg.enabled:
            # No-op; preserve legacy behavior
            continue

        present_values, has_missing, wc_list = _collect_present_values(comp)

        normalized: dict[str, list[str]] = {}
        for ent in wc_list:
            vals = present_values.get(ent, set())
            if not vals:
                # Entirely absent: drop this entity
                continue
            if has_missing.get(ent, False) and _in_scope(ent, s_cfg.scope):
                vals = set(vals) | {s_cfg.label}
            normalized[ent] = sorted(vals)

        _set_component_entities(comp, normalized)

        # annotate component with snakenull rendering preferences if helpful later
        if hasattr(comp, "__dict__"):
            comp.snakenull_label = s_cfg.label
            comp.snakenull_include_prefix = s_cfg.include_prefix
        if isinstance(comp, MutableMapping):
            comp["snakenull_label"] = s_cfg.label
            comp["snakenull_include_prefix"] = s_cfg.include_prefix

    return inputs
