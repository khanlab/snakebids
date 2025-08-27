# Fix Test Logic Bug in Participant Filtering Test

## Summary

This PR fixes a bug in the test `test_participant_label_doesnt_filter_comps_when_subject_has_filter_no_wcard` where the test was creating a modified config but not using it in the actual test call, causing intermittent test failures.

## Problem

The test was experiencing flaky behavior - sometimes passing, sometimes failing. Investigation revealed that the test was:

1. Creating a `config` from the dataset using `create_snakebids_config(rooted)`
2. Modifying this config by adding subject filters to each component
3. But then calling `generate_inputs()` with the original `create_snakebids_config(rooted)` instead of the modified `config`

This meant the subject filters were never actually applied during testing, making the test's behavior non-deterministic and dependent on the randomly generated test data from Hypothesis.

## Solution

Fixed the test to use the modified `config` variable in the `generate_inputs()` call instead of recreating the config. This ensures that the subject filters are properly applied during the test, making the test behavior consistent and deterministic.

## Changes

- Modified `test_participant_label_doesnt_filter_comps_when_subject_has_filter_no_wcard` in `snakebids/tests/test_generate_inputs.py`
- Changed the `generate_inputs()` call to use the modified `config` instead of `create_snakebids_config(rooted)`

## Testing

- The previously flaky test now passes consistently
- All existing tests continue to pass
- Pre-commit hooks (ruff, pyright, codespell) all pass

## Impact

This fix ensures reliable test execution and eliminates a source of CI flakiness. The test now properly validates the intended behavior: that participant label filtering doesn't affect components when subject entities have explicit filters applied.
