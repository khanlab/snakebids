"""Tools to generate a Snakemake-based BIDS app."""

import json
import os
import pathlib
import subprocess
import argparse
import logging
import sys
import yaml
import itertools as it
import shutil

import snakemake
from snakemake.io import load_configfile
from colorama import Fore

from bids import config as bidsconfig
from snakebids.exceptions import ConfigError
from snakebids.output import prepare_output, retrofit_output, write_config_file, write_output_mode

# We define Path here in addition to pathlib to put both variables in globals()
# This way, users specifying a path type in their config.yaml can indicate
# either Path or pathlib.Path
Path = pathlib.Path

bidsconfig.set_option("extension_initial_dot", True)
logger = logging.Logger(__name__)




class KeyValue(argparse.Action):
    """Class for accepting key=value pairs in argparse"""

    # Constructor calling
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())

        for value in values:
            # split it into key and value
            key, value = value.split("=")
            # assign into dictionary
            getattr(namespace, self.dest)[key] = value


class SnakemakeHelpAction(argparse.Action):
    """Class for printing snakemake usage in argparse"""

    def __call__(self, parser, namespace, values, option_string=None):
        run("snakemake -h")
        sys.exit(0)


def run(command, env=None):
    """Helper function for running a system command while merging
    stderr/stdout to stdout.

    Parameters
    ----------
    command : list of str
        command to run

    env : dict, optional
        environment variable to set before running the command
    """
    if env is None:
        env = {}

    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        env=merged_env,
    )
    while True:
        line = process.stdout.readline()
        line = str(line, "utf-8")[:-1]
        print(line)
        if line == "" and process.poll() is not None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d" % process.returncode)



def resolve_path(path_candidate):
    """Helper function to resolve any paths or list
    of paths it's passed. Otherwise, returns the argument 
    unchanged. 

    Parameters
    ----------
    command : list, os.Pathlike, object
        command to run

    Returns
    -------
    list, os.Pathlike, object
        If os.Pathlike or list  of os.Pathlike, the same paths resolved.
        Otherwise, the argument unchanged.
    """
    if isinstance(path_candidate, list):
        return [resolve_path(p) for p in path_candidate]

    if isinstance(path_candidate, os.PathLike):
        return path_candidate.resolve()
    
    return path_candidate

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
    "config/config.yaml"
]


class SnakeBidsApp:
    """Snakebids app with config and arguments.

    Parameters
    ----------
    snakemake_dir : str
        Root directory of the snakebids app, containing the config file and
        workflow files.
    skip_parse_args : bool, optional
        If true, the Snakebids app will not attempt to parse input arguments,
        and will only handle the config file.
    out_configfile : str
        Path to the updated configfile (YAML or JSON), relative to the output 
        working directory. This should be the same as the `configfile: ` used
        in your workflow. (default: 'config/snakebids.yml')

    Attributes
    ----------
    config : dict
        Contains all the configuration variables parsed from the config file
        and generated during the initialization of the SnakeBidsApp.
    parser_include_snakemake : ArgumentParser
        Parser including the generic Snakemake parser as a parent. This will
        contain all arguments a Snakemake app can receive.
    parser : ArgumentParser
        Parser including only the arguments specific to this Snakebids app, as
        specified in the config file.
    snakefile : str
        Absolute path to the input Snakefile
            join(snakemake_dir, snakefile_path)
    configfile_path : str
        Relative path to config file (relative to snakemake_dir)
    updated_config : str
        Absolute path to the updated config file to write.
    workflow_mode : bool
        If true, produces output in workflow mode, otherwise snakemake mode.
    force : bool
        If false, prevents conversion of output from workflow mode to bidsapp mode.
    outputdir : Path
        Path to outputdir specified by user
    retrofit : bool
        If true, run_snakemake will attempt to convert legacy output format into bidsapp
        format.
    """

    outputdir: Path
    snakemake_dir: Path

    def __init__(self, snakemake_dir, skip_parse_args=False):
        # input argument is the dir where snakemake would be run
        # we use this to locate the config file, and snakefile adding them to
        # generated_config and also add the snakemake_dir to the
        # generated_config, so the workflow can use it to source files from
        # it (e.g. atlases etc..)

        self.snakemake_dir = Path(snakemake_dir).resolve()

        # look for snakebids.yml in the snakemake_dir, quit if not found
        self.configfile_path = ""
        for path in CONFIGFILE_CHOICES:
            if (self.snakemake_dir / path).exists():
                self.configfile_path = path
                break
        if not self.configfile_path:
            raise ConfigError(
                f"Error: no config file found, tried {', '.join(CONFIGFILE_CHOICES)}."
            )

        # look for snakefile in the snakemake_dir, quit if not found
        self.snakefile = None
        for snakefile_path in SNAKEFILE_CHOICES:
            if (self.snakemake_dir / snakefile_path).exists():
                self.snakefile = Path(snakemake_dir, snakefile_path)
                break
        if self.snakefile is None:
            raise ConfigError(
                f"Error: no Snakefile found, tried {', '.join(SNAKEFILE_CHOICES)}."
            )

        self.config = load_configfile(Path(snakemake_dir,
                                            self.configfile_path))

        if self.config.get("debug", False):
            logging.basicConfig(level=logging.DEBUG)

        # add path to snakefile to the config -- so workflows can grab files
        # relative to the snakefile folder
        self.config["snakemake_dir"] = snakemake_dir
        self.config["snakefile"] = self.snakefile

        self.parser_include_snakemake = self._create_parser(
            include_snakemake=True
        )
        self.parser = self._create_parser()

        if not skip_parse_args:
            self._parse_args()


    def _create_parser(self, include_snakemake=False):
        """Create a parser with snakemake parser as parent solely for
        displaying help and checking conflicts, but then for actual parsing
        use snakebids parser to parse known args, then pass remaining to
        snakemake.
        """

        if include_snakemake:
            # get snakemake parser
            smk_parser = snakemake.get_argument_parser()

            # create parser
            parser = argparse.ArgumentParser(
                description="Snakebids helps build BIDS Apps with Snakemake",
                add_help=False,
                parents=[smk_parser],
            )
        else:
            parser = argparse.ArgumentParser(
                description="Snakebids helps build BIDS Apps with Snakemake"
            )

        parser.add_argument(
            "--workflow-mode",
            "-W",
            action="store_true"
        )

        # We use -x as the alias because both -f and -F are taken by snakemake
        parser.add_argument(
            "--force-conversion",
            "-x",
            action="store_true"
        )

        parser.add_argument(
            "--retrofit",
            action="store_true"
        )

        # add option for printing out snakemake usage
        parser.add_argument(
            "--help_snakemake",
            nargs=0,
            action=SnakemakeHelpAction,
            help=(
                "Options to Snakemake can also be passed directly at the "
                "command-line, use this to print Snakemake usage"
            ),
        )

        # create parser group for app options
        app_group = parser.add_argument_group(
            "SNAKEBIDS", "Options for snakebids app"
        )

        # update the parser with config options
        for name, parse_args in self.config["parse_args"].items():
            # Convert type annotations from strings to class types
            # We first check that the type annotation is, in fact,
            # a str to allow the edge case where it's already
            # been converted
            if "type" in parse_args and isinstance(parse_args["type"], str):
                try:
                    parse_args["type"] = globals()[parse_args["type"]]
                except KeyError as err:
                    raise TypeError(
                        f"{parse_args['type']} is not available "
                        + f"as a type for {name}"
                    ) from err
                    
            app_group.add_argument(name, **parse_args)

        # general parser for
        # --filter_{input_type} {key1}={value1} {key2}={value2}...
        # create filter parsers, one for each input_type
        filter_opts = parser.add_argument_group(
            "BIDS FILTERS",
            "Filters to customize PyBIDS get() as key=value pairs",
        )

        for input_type in self.config["pybids_inputs"].keys():
            argname = f"--filter_{input_type}"
            arglist_default = [
                f"{key}={value}"
                for (key, value) in self.config["pybids_inputs"][input_type][
                    "filters"
                ].items()
            ]
            arglist_default_string = " ".join(arglist_default)

            filter_opts.add_argument(
                argname,
                nargs="+",
                action=KeyValue,
                help=f"(default: {arglist_default_string})",
            )

        # general parser for
        # --wildcards_{input_type} {wildcard1} {wildcard2} ...
        # create wildcards parsers, one for each input_type
        wildcards_opts = parser.add_argument_group(
            "INPUT WILDCARDS",
            "File path entities to use as wildcards in snakemake",
        )

        for input_type in self.config["pybids_inputs"].keys():
            argname = f"--wildcards_{input_type}"
            arglist_default = [
                f"{wc}"
                for wc in self.config["pybids_inputs"][input_type]["wildcards"]
            ]
            arglist_default_string = " ".join(arglist_default)

            wildcards_opts.add_argument(
                argname,
                nargs="+",
                help=f"(default: {arglist_default_string})",
            )

        override_opts = parser.add_argument_group(
            "PATH OVERRIDE",
            (
                "Options for overriding BIDS by specifying absolute paths "
                "that include wildcards, e.g.: "
                "/path/to/my_data/{subject}/t1.nii.gz"
            ),
        )

        # create path override parser
        for input_type in self.config["pybids_inputs"].keys():
            argname = f"--path_{input_type}"
            override_opts.add_argument(argname, default=None)

        return parser

    def _parse_args(self):

        # use snakebids parser to parse the known arguments
        # will pass the rest of args when running snakemake
        all_args = self.parser.parse_known_args()

        args = all_args[0]
        snakemake_args = all_args[1]
        # resolve all path items to get absolute paths
        args.__dict__ = {
            k: resolve_path(v) for k, v in args.__dict__.items()
        }

        self.workflow_mode = args.workflow_mode
        self.force = args.force_conversion
        self.outputdir = Path(args.output_dir).resolve()
        self.retrofit = args.retrofit

        # add snakebids arguments to config
        self.config.update(args.__dict__)

        # add snakemake arguments to config
        self.config.update({"snakemake_args": snakemake_args})

        # argparse adds filter_{input_type} to the config
        # we want to update the pybids_inputs dict with this, then remove the
        # filter_{input_type} dict
        for input_type in self.config["pybids_inputs"].keys():
            arg_filter_dict = self.config[f"filter_{input_type}"]
            if arg_filter_dict is not None:
                self.config["pybids_inputs"][input_type]["filters"].update(
                    arg_filter_dict
                )
            del self.config[f"filter_{input_type}"]

        # add cmdline defined wildcards from the list:
        # wildcards_{input_type}
        for input_type in self.config["pybids_inputs"].keys():
            wildcards_list = self.config[f"wildcards_{input_type}"]
            if wildcards_list is not None:
                self.config["pybids_inputs"][input_type][
                    "wildcards"
                ] += wildcards_list
            del self.config[f"wildcards_{input_type}"]

        # add custom input paths to
        # config['pybids_inputs'][input_type]['custom_path']
        for input_type in self.config["pybids_inputs"].keys():
            custom_path = self.config[f"path_{input_type}"]
            if custom_path is not None:
                self.config["pybids_inputs"][input_type][
                    "custom_path"
                ] = Path(custom_path).resolve()
            del self.config[f"path_{input_type}"]

        # replace paths with realpaths
        self.config["bids_dir"] = Path(self.config["bids_dir"]).resolve()
        self.config["output_dir"] = Path(self.config["output_dir"]).resolve()

    
    def run_snakemake(self):
        """Run snakemake with that config.

        Workflow snakefile will read snakebids config, create inputs_config,
        and read that in.
        """

        if self.snakemake_dir == self.outputdir:
            write_output_mode(self.outputdir/".snakebids", "workflow")
            mode = "workflow"
            force_config_overwrite = True
            if self.workflow_mode == False:
                print(
                    f"{Fore.YELLOW}You specified your output to be in the snakebids "
                    "directory, so we're switching automatically to workflow mode!\n"

                    f"{Fore.RESET}You'll find your results in the "
                    "`snakebidsdir/results` directory."
                )
        else:
            mode = "workflow" if self.workflow_mode else "bidsapp"
            force_config_overwrite = False
            if self.retrofit:
                for path in CONFIGFILE_CHOICES:
                    if (self.outputdir / path).exists():
                        self.configfile_path = path
                if not retrofit_output(
                    self.outputdir,
                    # Find all config files in the outputdir
                    (
                        self.outputdir/p 
                        for p in CONFIGFILE_CHOICES if (self.outputdir/p).exists()
                    )
                ):
                    exit(1)
        
        root = prepare_output(
            self.snakemake_dir,
            self.outputdir,
            mode,
            self.force
        )

        if mode == "workflow":
            self.config["output_dir"] = str(root)
            new_config_file = self.outputdir/self.configfile_path
            
            cwd = self.outputdir
        else:
            new_config_file = self.outputdir/"code"/Path(self.configfile_path).name
            cwd = self.snakemake_dir
            force_config_overwrite = True

        write_config_file(new_config_file, self.config, force_config_overwrite)

        # running the chosen participant level
        analysis_level = self.config["analysis_level"]

        # run snakemake (passing any leftover args from argparse)
        # Filter any blank strings before submitting
        snakemake_argv = [*filter(None, [
            "--snakefile",
            str(self.snakefile),
            "--directory",
            str(cwd),
            "--configfile",
            str(new_config_file.resolve()),
            *self.config["snakemake_args"],
            *self.config["targets_by_analysis_level"][analysis_level],
        ])]

        snakemake.main(snakemake_argv)
        
