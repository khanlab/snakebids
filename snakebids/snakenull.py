from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, MutableMapping, Any

"""
snakenull.py: Post-processing utilities for normalizing mixed/absent entities
in Snakebids inputs using a configurable placeholder label (default: "snakenull").

This module operates on public objects returned by snakebids.generate_inputs()
and does not require internal hooks. You can later integrate it deeper if desired.
"""

@dataclass(frozen=True)
class SnakenullConfig:
    enabled: bool = False
    label: str = "snakenull"
    include_prefix: bool = True   # downstream path builders may use this
    scope: Iterable[str] | str = "all"  # "all" or iterable of entity names

def _in_scope(ent: str, scope: Iterable[str] | str) -> bool:
    return scope == "all" or ent in set(scope)

def _merge_cfg(global_cfg: Mapping[str, Any] | None,
               local_cfg: Mapping[str, Any] | None) -> SnakenullConfig:
    base = {"enabled": False, "label": "snakenull", "include_prefix": True, "scope": "all"}
    if global_cfg:
        base.update(global_cfg)
    if local_cfg:
        base.update(local_cfg)
    return SnakenullConfig(**base)

def _collect_present_values(component) -> tuple[Dict[str, set[str]], Dict[str, bool], list[str]]:
    """
    Work against both attribute- and mapping-style components.
    We infer:
      - wildcard list
      - matched files
      - each file's .entities mapping (PyBIDS-style)
    """
    # Wildcards list
    wc_list: list[str] = []
    if hasattr(component, "requested_wildcards"):
        wc_list = list(getattr(component, "requested_wildcards"))
    elif hasattr(component, "wildcards") and isinstance(component.wildcards, Mapping):
        wc_list = list(component.wildcards.keys())
    elif isinstance(component, Mapping) and "wildcards" in component:
        wc_list = list(component["wildcards"].keys())

    present_values: Dict[str, set[str]] = {w: set() for w in wc_list}
    has_missing: Dict[str, bool] = {w: False for w in wc_list}

    # Matched files
    files = None
    for attr in ("matched_files", "files", "items", "records"):
        if hasattr(component, attr):
            files = getattr(component, attr)
            break
    if files is None and isinstance(component, Mapping):
        for key in ("matched_files", "files", "items", "records"):
            if key in component:
                files = component[key]
                break
    files = files or []

    for rec in files:
        ents = {}
        if hasattr(rec, "entities"):
            ents = getattr(rec, "entities") or {}
        elif isinstance(rec, Mapping) and "entities" in rec:
            ents = rec["entities"] or {}
        for ent in wc_list:
            if ent in ents and ents[ent] not in (None, ""):
                present_values[ent].add(str(ents[ent]))
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
    wildcards = {k: "{" + k + "}" for k in entities.keys()}
    if hasattr(component, "wildcards"):
        component.wildcards = wildcards
    elif isinstance(component, MutableMapping):
        component["wildcards"] = wildcards
    # optional flags for downstream path builders
    if hasattr(component, "__dict__"):
        setattr(component, "_snakenull_label", None)
        setattr(component, "_snakenull_include_prefix", True)
    if isinstance(component, MutableMapping):
        component["_snakenull_label"] = None
        component["_snakenull_include_prefix"] = True

def normalize_inputs_with_snakenull(
    inputs: Mapping[str, Any],
    *,
    config: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    """
    Post-process the result of snakebids.generate_inputs() to:
      1) Skip entities listed in wildcards that are entirely absent in the dataset
      2) For entities present in some files but missing in others, insert a placeholder
         value (default: "snakenull")

    Parameters
    ----------
    inputs : Mapping[str, Any]
        The object returned by snakebids.generate_inputs(), mapping
        component name -> component.
    config : Mapping[str, Any], optional
        Top-level Snakebids config dict. Read from:
           config["snakenull"] (global defaults)
           config["pybids_inputs"][<component>]["snakenull"] (per-component override)

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

        normalized: Dict[str, list[str]] = {}
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
            setattr(comp, "_snakenull_label", s_cfg.label)
            setattr(comp, "_snakenull_include_prefix", s_cfg.include_prefix)
        if isinstance(comp, MutableMapping):
            comp["_snakenull_label"] = s_cfg.label
            comp["_snakenull_include_prefix"] = s_cfg.include_prefix

    return inputs
