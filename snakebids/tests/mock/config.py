# ruff: noqa

from __future__ import annotations

from snakebids.types import InputsConfig

pybids_inputs: InputsConfig = {
    "bold": {
        "filters": {"suffix": "bold", "extension": ".nii.gz", "datatype": "func"},
        "wildcards": ["subject", "session", "acquisition", "task", "run"],
    }
}

parse_args = {
    "bids_dir": {
        "help": "The directory with the input dataset formatted according\nto the BIDS standard."
    },
    "output_dir": {
        "help": "The directory where the output files\nshould be stored. If you are running group level analysis\nthis folder should be prepopulated with the results of the\nparticipant level analysis."
    },
    "analysis_level": {
        "help": "Level of the analysis that will be performed.",
        "choices": ["participant"],
    },
    "--participant_label": {
        "help": 'The label(s) of the participant(s) that should be analyzed. The label\ncorresponds to sub-<participant_label> from the BIDS spec \n(so it does not include "sub-"). If this parameter is not \nprovided all subjects should be analyzed. Multiple \nparticipants can be specified with a space separated list.',
        "nargs": "+",
    },
    "--exclude_participant_label": {
        "help": 'The label(s) of the participant(s) that should be excluded. The label\ncorresponds to sub-<participant_label> from the BIDS spec \n(so it does not include "sub-"). If this parameter is not \nprovided all subjects should be analyzed. Multiple \nparticipants can be specified with a space separated list.',
        "nargs": "+",
    },
    "--derivatives": {
        "help": "Path(s) to a derivatives dataset, for folder(s) that contains multiple derivatives datasets (default: %(default)s) ",
        "default": False,
        "type": "Path",
        "nargs": "+",
    },
    "--skip_bids_validation": {
        "help": "Skip validation of BIDS dataset. BIDS validation is performed by "
        "default using the bids-validator (if installed) or with the pybids "
        "validator implementation (if bids-validator is not installed). "
        "(default: %(default)s)",
        "dest": "plugins.validator.skip",
        "action": "store_true",
        "default": False,
    },
    "--arg-using-dash-syntax": {
        "help": "A fake argument for testing purposes",
        "nargs": "+",
    },
}

config = {
    "pybids_inputs": pybids_inputs,
    "targets_by_analysis_level": {"participant": [""]},
    "analysis_levels": ["participant"],
    "parse_args": parse_args,
}
