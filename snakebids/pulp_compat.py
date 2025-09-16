"""
Compatibility layer for PuLP API changes.

This module ensures compatibility between different versions of PuLP that may use
different method names for listing solvers (listSolvers vs list_solvers).
"""

from __future__ import annotations

import warnings


def ensure_pulp_compatibility() -> None:
    """
    Ensure compatibility between different PuLP versions.

    This function patches the pulp module to provide both listSolvers and list_solvers
    methods, regardless of which version is installed.
    """
    try:
        import pulp  # type: ignore[import-untyped]

        # Check if list_solvers exists but listSolvers doesn't (newer PuLP)
        if hasattr(pulp, "list_solvers") and not hasattr(pulp, "listSolvers"):
            # Add the camelCase version for backward compatibility
            pulp.listSolvers = pulp.list_solvers  # type: ignore[attr-defined]

        # Check if listSolvers exists but list_solvers doesn't (older PuLP)
        elif hasattr(pulp, "listSolvers") and not hasattr(pulp, "list_solvers"):
            # Add the snake_case version for forward compatibility
            pulp.list_solvers = pulp.listSolvers  # type: ignore[attr-defined]

        # If neither exists, provide a fallback
        elif not hasattr(pulp, "listSolvers") and not hasattr(pulp, "list_solvers"):

            def _fallback_list_solvers(
                only_available: bool = True,
            ) -> list[str]:
                """Fallback solver list if PuLP doesn't provide one."""
                warnings.warn(
                    "PuLP solver listing not available. Using fallback list.",
                    UserWarning,
                    stacklevel=2,
                )
                return ["PULP_CBC_CMD"]

            pulp.listSolvers = _fallback_list_solvers  # type: ignore[attr-defined]
            pulp.list_solvers = _fallback_list_solvers  # type: ignore[attr-defined]

    except (ImportError, AttributeError, TypeError, RuntimeError) as e:
        # Handle expected errors during PuLP compatibility patching
        warnings.warn(
            f"Failed to apply PuLP compatibility patch: {e}. "
            "This may cause issues with ILP scheduling.",
            UserWarning,
            stacklevel=2,
        )


# Auto-apply compatibility when module is imported
ensure_pulp_compatibility()
