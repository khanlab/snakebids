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
    scope:
        Either ``"all"`` or an iterable of entity names to normalize.
    """

    enabled: bool = False
    label: str = "snakenull"
    scope: Iterable[str] | str = "all"  # "all" or iterable of entity names


def _in_scope(ent: str, scope: Iterable[str] | str) -> bool:
    """Return True if an entity name is within the configured normalization scope."""
    return scope == "all" or ent in set(scope)


def _merge_cfg(
    global_cfg: Mapping[str, Any] | None, local_cfg: Mapping[str, Any] | None
) -> SnakenullConfig:
    """Merge global and per-component snakenull configuration into a dataclass."""
    base: dict[str, Any] = {
        "enabled": False,
        "label": "snakenull",
        "scope": "all",
    }
    if global_cfg:
        base.update(global_cfg)
    if local_cfg:
        base.update(local_cfg)

    # Type-safe extraction with defaults
    enabled = bool(base.get("enabled", False))
    label = str(base.get("label", "snakenull"))
    scope = base.get("scope", "all")

    return SnakenullConfig(enabled=enabled, label=label, scope=scope)


def _wildcards_list_from_component(component: Any) -> list[str]:
    """Return wildcard names from the component in a tolerant way."""

    def safe_convert_to_strings(items: Any) -> list[str]:
        try:
            return [str(item) for item in items]
        except (TypeError, ValueError):
            return []

    # Try requested_wildcards attribute first
    if hasattr(component, "requested_wildcards"):
        requested = getattr(component, "requested_wildcards", None)
        if requested is not None:
            result = safe_convert_to_strings(requested)
            if result:
                return result

    # Try wildcards attribute
    if hasattr(component, "wildcards"):
        wildcards = getattr(component, "wildcards", None)
        if isinstance(wildcards, Mapping):
            result = safe_convert_to_strings(wildcards.keys())
            if result:
                return result

    # Try wildcards in mapping-like component
    if isinstance(component, Mapping) and "wildcards" in component:
        w: Any = component["wildcards"]  # type: ignore[misc]
        if isinstance(w, Mapping):
            result = safe_convert_to_strings(w.keys())
            if result:
                return result
        if isinstance(w, (list, tuple, set)):
            result = safe_convert_to_strings(w)
            if result:
                return result

    return []


# Add above `_collect_present_values`
def _collect_from_records(
    records: list[Mapping[str, Any]], wc_list: list[str]
) -> tuple[dict[str, set[str]], dict[str, bool]]:
    """Accumulate present/missing from a list of per-file entity dicts."""
    present_values: dict[str, set[str]] = {w: set() for w in wc_list}
    has_missing: dict[str, bool] = {w: False for w in wc_list}
    for rec in records:
        get = rec.get
        for ent in wc_list:
            val = get(ent)
            if val not in (None, ""):
                present_values[ent].add(str(val))
            else:
                has_missing[ent] = True
    return present_values, has_missing


def _collect_from_zip_lists(
    zip_lists: Mapping[str, list[str]], wc_list: list[str]
) -> tuple[dict[str, set[str]], dict[str, bool]]:
    """Accumulate present/missing from rectangular zip_lists."""
    present_values: dict[str, set[str]] = {w: set() for w in wc_list}
    has_missing: dict[str, bool] = {w: False for w in wc_list}
    total = max((len(zip_lists.get(ent, [])) for ent in wc_list), default=0)
    for ent in wc_list:
        lst = list(zip_lists.get(ent, []))
        if len(lst) < total:
            lst.extend([""] * (total - len(lst)))
        for val in lst:
            if val not in (None, ""):
                present_values[ent].add(str(val))
            else:
                has_missing[ent] = True
    return present_values, has_missing


def _collect_present_values(
    component: Any,
) -> tuple[dict[str, set[str]], dict[str, bool], list[str]]:
    """Collect present values and detect missing ones."""
    wc_list: list[str] = _wildcards_list_from_component(component)

    # 1) Primary: precomputed per-file dicts (added in _get_component)
    records = getattr(component, "snakenull_entities_per_file", None)
    if records:
        present, missing = _collect_from_records(list(records), wc_list)
        return present, missing, wc_list

    # 2) Fallback: matched_files with entities attribute (test components)
    if hasattr(component, "matched_files"):
        matched_files = getattr(component, "matched_files", None)
        if matched_files:
            records_from_files: list[Mapping[str, Any]] = []
            for file_obj in matched_files:
                if hasattr(file_obj, "entities"):
                    entities = getattr(file_obj, "entities", {})
                    if isinstance(entities, Mapping):
                        records_from_files.append(entities)  # type: ignore[arg-type]
            if records_from_files:
                present, missing = _collect_from_records(records_from_files, wc_list)
                return present, missing, wc_list

    # 3) Fallback: zip_lists (e.g., custom_path components)
    zip_lists: Mapping[str, list[str]] | None = None
    if hasattr(component, "zip_lists"):
        zip_lists_attr = getattr(component, "zip_lists", None)
        if isinstance(zip_lists_attr, Mapping):
            zip_lists = zip_lists_attr  # type: ignore[assignment]
    elif isinstance(component, Mapping):
        zl: Any = component.get("zip_lists")  # type: ignore[misc]
        if isinstance(zl, Mapping):
            zip_lists = zl  # type: ignore[assignment]
    if zip_lists:
        present, missing = _collect_from_zip_lists(zip_lists, wc_list)
        return present, missing, wc_list

    # 4) No data available
    present_values: dict[str, set[str]] = {w: set() for w in wc_list}
    has_missing: dict[str, bool] = {w: False for w in wc_list}
    return present_values, has_missing, wc_list


def _set_component_entities(component: Any, entities: Mapping[str, list[str]]) -> None:
    """Expose normalized entity domains without breaking frozen components.

    - If the component is a MutableMapping, write into its keys
    ('entities', 'wildcards').
    - Otherwise, *best effort* set attributes, but swallow AttributeError
    from frozen attrs.
    - Always try to stash under a side attribute 'snakenull_entities' for
    consumers that
      choose to read it explicitly.
    """
    # Prefer writing into mapping-like components (safe, no __setattr__)
    if isinstance(component, MutableMapping):
        component["entities"] = dict(entities)
        component["wildcards"] = {k: "{" + k + "}" for k in entities}
        # Optional hints for downstream renderers
        existing_label = None
        existing_prefix = True
        try:
            # Use Any cast to avoid type checker issues with dynamic getattr
            comp_any: Any = component  # type: ignore[misc]
            existing_label = getattr(comp_any, "_snakenull_label", None)
            existing_prefix = getattr(comp_any, "_snakenull_include_prefix", True)
        except (AttributeError, TypeError):
            pass
        component["_snakenull_label"] = existing_label
        component["_snakenull_include_prefix"] = existing_prefix
        return

    # Attribute-based components: be defensive (attrs-frozen will raise on __setattr__)
    import contextlib

    with contextlib.suppress(Exception):
        if hasattr(component, "entities"):
            component.entities = dict(entities)
    with contextlib.suppress(Exception):
        if hasattr(component, "wildcards"):
            component.wildcards = {k: "{" + k + "}" for k in entities}

    # Side channel (optional): stash normalized domains for readers that opt-in
    with contextlib.suppress(Exception):
        component.snakenull_entities = dict(entities)


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
    pybids_inputs: Mapping[str, Any] = {}
    if config and isinstance(config.get("pybids_inputs"), Mapping):
        pybids_inputs = config["pybids_inputs"]
    global_cfg: Mapping[str, Any] = {}
    if config and isinstance(config.get("snakenull"), Mapping):
        global_cfg = config["snakenull"]

    for cname, comp in inputs.items():
        local_cfg: Mapping[str, Any] = {}
        comp_cfg = pybids_inputs.get(cname)
        if isinstance(comp_cfg, Mapping):
            # Use Any to handle unknown mapping types safely
            comp_cfg_any: Any = comp_cfg  # type: ignore[misc]
            snakenull_cfg: Any = comp_cfg_any.get("snakenull", {})
            if snakenull_cfg and isinstance(snakenull_cfg, Mapping):
                local_cfg = snakenull_cfg  # type: ignore[assignment]
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
            comp.snakenull_label = s_cfg.label  # type: ignore[attr-defined]
        if isinstance(comp, MutableMapping):
            comp["snakenull_label"] = s_cfg.label

    return inputs
