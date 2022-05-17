"""Tools to generate a Snakemake-based BIDS app."""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import attr
import boutiques.creator as bc
import snakemake
from snakemake.io import load_configfile

from snakebids.cli import (
    SnakebidsArgs,
    add_dynamic_args,
    create_parser,
    parse_snakebids_args,
)
from snakebids.exceptions import ConfigError, RunError
from snakebids.utils.output import (
    prepare_bidsapp_output,
    write_config_file,
    write_output_mode,
)

logger = logging.Logger(__name__)


SNAKEFILE_CHOICES = [
    "Snakefile",
    "snakefile",
    "workflow/Snakefile",
    "workflow/snakefile",
]


CONFIGFILE_CHOICES = [
    "config/snakebids.yml",
    "config/snakebids.yaml",
    "config/snakebids.json",
    "snakebids.yml",
    "snakebids.yaml",
    "snakebids.json",
    "config.yml",
    "config.yaml",
    "config.json",
    "config/config.json",
    "config/config.yml",
    "config/config.yaml",
]


def _get_file_paths(choices: List[str], file_name: str):
    def wrapper(self: "SnakeBidsApp"):
        for path in choices:
            if (self.snakemake_dir / path).exists():
                if file_name == "config":
                    return Path(path)
                # else, snakefile
                return Path(self.snakemake_dir, path)

        raise ConfigError(
            f"Error: no {file_name} file found, tried {', '.join(choices)}."
        )

    return wrapper


# pylint: disable=unsubscriptable-object, unsupported-assignment-operation,
# pylint: disable=too-few-public-methods
@attr.define(slots=False)
class SnakeBidsApp:
    """Snakebids app with config and arguments.

    Attributes
    ----------
    snakemake_dir : str
        Root directory of the snakebids app, containing the config file and workflow
        files.
    parser : ArgumentParser, optional
        Parser including only the arguments specific to this Snakebids app, as specified
        in the config file. By default, it will use `create_parser()` from `cli.py`
    configfile_path : str, optional
        Relative path to config file (relative to snakemake_dir). By default,
        autocalculates based on snamake_dir
    snakefile_path : str, optional
        Absolute path to the input Snakefile. By default, autocalculates based on
        snakemake_dir::

            join(snakemake_dir, snakefile_path)
    config : dict, optional
        Contains all the configuration variables parsed from the config file
        and generated during the initialization of the SnakeBidsApp.
    args : SnakebidsArgs, optional
        Arguments to use when running the app. By default, generated using the parser
        attribute, autopopulated with args from `config.py`
    """

    snakemake_dir: Path = attr.ib(converter=lambda path: Path(path).resolve())
    skip_parse_args: bool = False
    parser: argparse.ArgumentParser = create_parser()
    configfile_path: Path = attr.Factory(
        _get_file_paths(CONFIGFILE_CHOICES, "config"), takes_self=True
    )
    snakefile_path: Path = attr.Factory(
        _get_file_paths(SNAKEFILE_CHOICES, "Snakefile"), takes_self=True
    )
    config: Dict[str, Any] = attr.Factory(
        lambda self: load_configfile(self.snakemake_dir / self.configfile_path),
        takes_self=True,
    )
    args: Optional[SnakebidsArgs] = None

    def run_snakemake(self):
        """Run snakemake with that config.

        Workflow snakefile will read snakebids config, create inputs_config,
        and read that in.
        """

        # If no SnakebidsArgs were provided on class instantiation, we compute args
        # using the provided parser
        if self.args:
            args = self.args
        else:
            # Dynamic args include --filter-... and --wildcards-... . They depend on the
            # config
            add_dynamic_args(
                self.parser, self.config["parse_args"], self.config["pybids_inputs"]
            )
            args = parse_snakebids_args(self.parser)

        # Update our config file:
        # - Add path to snakefile to the config so workflows can grab files relative to
        #    the snakefile folder
        # - Add info from args
        # - Set mode (bidsapp or workflow) and output_dir appropriately
        self.config["snakemake_dir"] = self.snakemake_dir
        self.config["snakefile"] = self.snakefile_path

        # Update config with pybids settings
        if args.pybidsdb_dir:
            self.config["pybids_db_dir"] = args.pybidsdb_dir
        self.config["pybids_db_reset"] = args.reset_db

        update_config(self.config, args)

        # First, handle outputs in snakebids_root or results folder
        try:
            # py3.9 has the Path.is_relative() function. But as long as we support py38
            # and lower, this is the easiest way
            args.outputdir.resolve().relative_to(self.snakemake_dir / "results")
            relative_to_results = True
        except ValueError:
            relative_to_results = False

        if self.snakemake_dir == args.outputdir.resolve() or relative_to_results:
            write_output_mode(self.snakemake_dir / ".snakebids", "workflow")

            new_config_file = self.snakemake_dir / self.configfile_path
            cwd = self.snakemake_dir

            if self.config["output_dir"] == self.snakemake_dir.resolve():
                self.config["output_dir"] /= "results"
                self.config["root"] = "results"
                # Print a friendly warning that the output directory will change
                logger.info(
                    "You specified your output to be in the snakebids directory, so "
                    "we're automatically putting your outputs in the results "
                    "subdirectory.\nYou'll find your results in `%s`",
                    (self.snakemake_dir / "results").resolve(),
                )
            else:
                self.config["root"] = ""

        # else, we run in bidsapp mode
        else:
            # Attempt to prepare the output folder. Anything going wrong will raise a
            # RunError, as described in the docstring
            try:
                prepare_bidsapp_output(args.outputdir, args.force)
            except RunError as err:
                print(err.msg)
                sys.exit(1)
            cwd = args.outputdir
            new_config_file = args.outputdir / self.configfile_path
            self.config["root"] = ""

        # Write the config file
        write_config_file(
            config_file=new_config_file,
            data=self.config,
            force_overwrite=True,
        )

        # Run snakemake (passing any leftover args from argparse)
        # Filter any blank strings before submitting
        snakemake.main(
            [
                *filter(
                    None,
                    [
                        "--snakefile",
                        str(self.snakefile_path),
                        "--directory",
                        str(cwd),
                        "--configfile",
                        str(new_config_file.resolve()),
                        *self.config["snakemake_args"],
                        *self.config["targets_by_analysis_level"][
                            self.config["analysis_level"]
                        ],
                    ],
                )
            ]
        )

    def create_descriptor(self, out_file):
        """Generate a boutiques descriptor for this Snakebids app."""
        new_descriptor = bc.CreateDescriptor(self.parser, execname="run.py")
        new_descriptor.save(out_file)


def update_config(config: Dict[str, Any], snakebids_args: SnakebidsArgs):
    # add snakemake arguments to config
    config.update({"snakemake_args": snakebids_args.snakemake_args})

    # argparse adds filter_{input_type}
    # we want to update the pybids_inputs dict with this, then remove the
    # filter_{input_type} dict
    pybids_inputs = config["pybids_inputs"]
    args = snakebids_args.args_dict
    for input_type in pybids_inputs.keys():
        arg_filter_dict = args[f"filter_{input_type}"]
        if arg_filter_dict is not None:
            pybids_inputs[input_type]["filters"].update(arg_filter_dict)
        del args[f"filter_{input_type}"]

    # add cmdline defined wildcards from the list:
    # wildcards_{input_type}
    for input_type in pybids_inputs.keys():
        wildcards_list = args[f"wildcards_{input_type}"]
        if wildcards_list is not None:
            pybids_inputs[input_type]["wildcards"] += wildcards_list
        del args[f"wildcards_{input_type}"]

    # add custom input paths to
    # config['pybids_inputs'][input_type]['custom_path']
    for input_type in pybids_inputs.keys():
        custom_path = args[f"path_{input_type}"]
        if custom_path is not None:
            pybids_inputs[input_type]["custom_path"] = Path(custom_path).resolve()
        del args[f"path_{input_type}"]

    # add snakebids arguments to config
    config.update(args)

    # replace paths with realpaths
    config["bids_dir"] = Path(config["bids_dir"]).resolve()
    config["output_dir"] = Path(config["output_dir"]).resolve()
