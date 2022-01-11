"""Tools to generate a Snakemake-based BIDS app."""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import attr
import boutiques.creator as bc
import snakemake
from colorama import Fore
from snakemake.io import load_configfile

from snakebids.cli import (
    SnakebidsArgs,
    add_dynamic_args,
    create_parser,
    parse_snakebids_args,
)
from snakebids.exceptions import ConfigError, RunError
from snakebids.utils.output import (
    prepare_output,
    retrofit_output,
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


def _get_file_paths(choices, file_name, snakemake_dir):
    for path in choices:
        if (snakemake_dir / path).exists():
            return Path(snakemake_dir, path)
    raise ConfigError(
        f"Error: no {file_name} file found, tried {', '.join(CONFIGFILE_CHOICES)}."
    )


# pylint: disable=unsubscriptable-object, unsupported-assignment-operation,
# pylint: disable=too-few-public-methods
@attr.define(slots=False)
class SnakeBidsApp:
    """Snakebids app with config and arguments.

    Attributes
    ----------
    snakemake_dir : Path
        Root directory of the snakebids app, containing the config file and workflow
        files.
    parser : ArgumentParser
        Parser including only the arguments specific to this Snakebids app, as specified
        in the config file. By default, it will use `create_parser()` from `cli.py`
    parser_run : ArgumentParser
        Subparser containing the run-specific arguments.
    configfile_path : str
        Relative path to config file (relative to snakemake_dir). By default,
        autocalculates based on snakemake_dir
    snakefile_path : str
        Absolute path to the input Snakefile. By default, autocalculates based on
        snakemake_dir
            join(snakemake_dir, snakefile_path)
    config : dict
        Contains all the configuration variables parsed from the config file
        and generated during the initialization of the SnakeBidsApp.
    args : SnakebidsArgs, optional
        Arguments to use when running the app. By default, generated using the parser
        attribute, autopopulated with args from `config.py`
    """

    snakemake_dir: Path = attr.ib(converter=lambda path: Path(path).resolve())
    parser: argparse.ArgumentParser
    parser_run: argparse.ArgumentParser
    configfile_path: Path
    snakefile_path: Path
    config: Dict[str, Any]
    args: Optional[SnakebidsArgs] = None
    skip_parse_args: bool = False

    @classmethod
    def from_filesystem(cls, snakemake_dir, args=None, skip_parse_args=False):
        """Build a SnakebidsApp by traversing a workflow on the filesystem.

        This will automatically locate the config file and Snakefile on the filesystem
        and construct the parser itself.

        Arguments
        ---------
        snakemake_dir : str
            Root directory of the snakebids app, containing the config file and workflow
            files.
        args : SnakebidsArgs, optional
            Arguments to use when running the app.
        skip_parse_args : bool
            The SnakebidsApp won't parse arguments if true.
        """
        parser, parser_run = create_parser()
        configfile_path = _get_file_paths(CONFIGFILE_CHOICES, "config", snakemake_dir)
        snakefile_path = _get_file_paths(CONFIGFILE_CHOICES, "Snakefile", snakemake_dir)
        config = load_configfile(snakemake_dir / configfile_path)
        return cls(
            snakemake_dir=snakemake_dir,
            parser=parser,
            parser_run=parser_run,
            configfile_path=configfile_path,
            snakefile_path=snakefile_path,
            config=config,
            args=args,
            skip_parse_args=skip_parse_args,
        )

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
                self.parser_run, self.config["parse_args"], self.config["pybids_inputs"]
            )
            args = parse_snakebids_args(self.parser)

        if hasattr(args, "path_boutiques"):
            self.create_descriptor(args.path_boutiques)
            print(f"Boutiques descriptor created at {args.path_boutiques}")
            return

        # If the snakemake_dir is the same as the outputdir, we need to switch into
        # workflow mode
        if self.snakemake_dir == args.outputdir:
            write_output_mode(args.outputdir / ".snakebids", "workflow")
            mode = "workflow"
            # The new config file will inevitably have the same path as the old, so we
            # allow overwriting
            force_config_overwrite = True

            # Print a friendly warning if the user didn't specify workflow mode
            if not args.workflow_mode:
                print(
                    f"{Fore.YELLOW}You specified your output to be in the snakebids "
                    "directory, so we're switching automatically to workflow mode!\n"
                    f"{Fore.RESET}You'll find your results in "
                    f"`{(self.snakemake_dir/'results').resolve()}`"
                )

        # Otherwise, both workflow and bidsapp mode are possible
        else:
            mode = "workflow" if args.workflow_mode else "bidsapp"
            # Disable config_overwrite to prevent accidental file modification
            # This will be set to true in the case of bidsapp mode
            force_config_overwrite = False

            # If the user asked to retrofit, we attempt it
            if args.retrofit and not retrofit_output(
                args.outputdir,
                # Find all config files in the outputdir
                (
                    args.outputdir / p
                    for p in CONFIGFILE_CHOICES
                    if (args.outputdir / p).exists()
                ),
            ):
                # If we get here, there was an error in retrofitting
                print(f"{Fore.YELLOW}Exiting. No conversion performed.{Fore.RESET}")
                sys.exit(1)

        # Attempt to prepare the output folder. Anything going wrong will raise a
        # RunError, as described in the docstring
        try:
            root = prepare_output(self.snakemake_dir, args.outputdir, mode, args.force)
        except RunError as err:
            print(err.msg)
            sys.exit(1)

        # Update our config file:
        # - Add path to snakefile to the config so workflows can grab files relative to
        #    the snakefile folder
        # - Add info from args
        # - Set mode (bidsapp or workflow) and output_dir appropriately
        self.config["snakemake_dir"] = self.snakemake_dir
        self.config["snakefile"] = self.snakefile_path
        update_config(self.config, args)
        if mode == "workflow":
            self.config["output_dir"] = str(root)
            self.config["root"] = "results"
        else:
            force_config_overwrite = True
            self.config["root"] = ""

        # Write the config file
        new_config_file = args.outputdir / self.configfile_path
        write_config_file(
            config_file=new_config_file,
            data=self.config,
            force_overwrite=force_config_overwrite,
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
                        str(args.outputdir),
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
