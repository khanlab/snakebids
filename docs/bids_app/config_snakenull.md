# Configuration: Snakenull normalization (optional)

Snakebids can **normalize mixed or absent entities**—cases where some files have an entity
like `acq` or `ses` but others don’t—by assigning a placeholder label to missing values
and skipping entities that are entirely absent. This stabilizes wildcard expansion and
output naming without requiring users to re-label inputs.

This behavior is **off by default** and can be controlled globally or per component.

## Global defaults

```yaml
snakenull:
  enabled: false       # default: legacy behavior
  label: snakenull     # placeholder value for missing entities
  scope: all           # "all" or a list, e.g., ["session", "acquisition"]
```
## Per-component override
```yaml
pybids_inputs:
  t1w:
    filters:
      suffix: T1w
      extension: .nii.gz
    wildcards: [subject, session, acquisition]
    snakenull:
      enabled: true
      scope: [session, acquisition]
```
## Usage
include these parameters in generate_inputs():
```python
inputs = generate_inputs(
    bids_dir=config["bids_dir"],
    pybids_inputs=config["pybids_inputs"],
    ...
    snakenull=config.get("snakenull"),
)
```
This mutates inputs in place so that:

- Entities listed in wildcards but entirely absent in your dataset are removed.

- Entities that are present for some files but missing for others are assigned the
placeholder label (default snakenull) and included in entities for expansion.

Downstream, your BidsComponent.entities will include e.g.:
```vbnet
acquisition: ["MPRAGE", "snakenull"]
session: ["01", "snakenull"]
```
so the expansion space (and target sets) remains stable even as new sessions or
acquisitions are added later.
