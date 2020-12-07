#!/usr/bin/env python3
import os
import subprocess
import sys
import yaml
from bids import BIDSLayout
import bids
from glob import glob
from snakemake import snakemake
import argparse
import json

bids.config.set_option('extension_initial_dot', True)

# create a keyvalue class for accepting key=value pairs in argparse
class keyvalue(argparse.Action): 
    # Constructor calling 
    def __call__( self , parser, namespace, 
                 values, option_string = None): 
        setattr(namespace, self.dest, dict()) 
          
        for value in values: 
            # split it into key and value 
            key, value = value.split('=') 
            # assign into dictionary 
            getattr(namespace, self.dest)[key] = value 

def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=True,
                               env=merged_env)
    while True:
        line = process.stdout.readline()
        line = str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)



SNAKEFILE_CHOICES = [
    "Snakefile",
    "snakefile",
    "workflow/Snakefile",
    "workflow/snakefile",
]

CONFIGFILE_CHOICES = [
    "config/snakebids.yml",
    "snakebids.yml"
]


class SnakeBidsApp:

    def __init__(self, snakemake_dir, skip_parse_args=False):

        #input argument is the dir where snakemake would be run
        # we use this to locate the config file, and snakefile adding them to generated_config
        # and also add the snakemake_dir to the generated_config, so the workflow can use it to 
        # source files from it (e.g. atlases etc..)
 
        #look for snakebids.yml in the snakemake_dir, quit if not found
        self.snakebids_config = None
        for p in CONFIGFILE_CHOICES:
            if os.path.exists(os.path.join(snakemake_dir,p)):
                self.snakebids_config = os.path.join(snakemake_dir,p)
                break
        if self.snakebids_config == None:
            print(
            "Error: no  snakebids.yml config file found, tried {}.".format(
                    ", ".join(SNAKEFILE_CHOICES), file=sys.stderr
                )
            )
            sys.exit(1)


        #look for snakefile in the snakemake_dir, quit if not found
        self.snakefile = None
        for p in SNAKEFILE_CHOICES:
            if os.path.exists(os.path.join(snakemake_dir,p)):
                self.snakefile = os.path.join(snakemake_dir,p)
                break
        if self.snakefile == None:
            print(
            "Error: no Snakefile found, tried {}.".format(
                    ", ".join(SNAKEFILE_CHOICES), file=sys.stderr
                )
            )
            sys.exit(1)



        self.__load_config()

        #add path to snakefile to the config -- so workflows can grab files relative to the snakefile folder  
        self.config['snakemake_dir'] = snakemake_dir

        self.parser = self.__create_parser()

        if not skip_parse_args:
            self.__parse_args()



    def __load_config(self):

        #load up workflow config file
        with open(self.snakebids_config, 'r') as infile:
            self.config = yaml.load(infile, Loader=yaml.FullLoader)




    def __create_parser(self):      

        #replace with proper name of pipeline here
        parser = argparse.ArgumentParser(description='snakebids-app')

        #update the parser with config options
        for name, parse_args in self.config['parse_args'].items():
            parser.add_argument(name, **parse_args)

        #general parser for --filter_{input_type} {key1}={value1} {key2}={value2}...
        #create filter parsers, one for each input_type

        for input_type in self.config['pybids_inputs'].keys():
            argname=f'--filter_{input_type}'
            arginstance=f'filter_{input_type}'.upper() #for help description
            arglist_default = [ f'{key}={value}'  for (key,value) in self.config['pybids_inputs'][input_type]['filters'].items() ]
            arglist_default_string = ' '.join(arglist_default)
 
            parser.add_argument(argname,nargs='+',action=keyvalue,
                                help=f'Filters (PyBIDS) for {input_type}, where {arginstance} is '
                                       f'key=value pair(s) (default: {arglist_default_string})')
 
        return parser


    def __parse_args(self):


        all_args = self.parser.parse_known_args()

        args = all_args[0]
        snakemake_args = all_args[1]


        #add arguments to config
        self.config.update(args.__dict__)
        self.config.update({'snakemake_args': snakemake_args})

        #argparse adds filter_{input_type} to the config
        # we want to update the pybids_inputs dict with this, then remove the filter_{input_type} dict
        for input_type in self.config['pybids_inputs'].keys():
            arg_filter_dict = self.config[f'filter_{input_type}']
            if arg_filter_dict is not None:
                self.config['pybids_inputs'][input_type]['filters'].update(arg_filter_dict)
            del self.config[f'filter_{input_type}']  
    

        #replace paths with realpaths
        self.config['bids_dir'] = os.path.realpath(self.config['bids_dir'])
        self.config['output_dir'] = os.path.realpath(self.config['output_dir'])
       


    def __write_updated_config(self):
        #create an updated snakebids config file
        self.updated_config = os.path.join(self.config['output_dir'],'config','snakebids.yml')

        #write it to file
        os.makedirs(os.path.dirname(self.updated_config),exist_ok=True)
        with open(self.updated_config, 'w') as outfile:
            yaml.dump(self.config, outfile, default_flow_style=False)




    # run snakemake with that config
    #workflow snakefile will read snakebids config, create inputs_config, read that in
    def run_snakemake(self):
        
        #write updated config
        self.__write_updated_config()

        # running the chosen participant level

        analysis_level = self.config['analysis_level']
        #runs snakemake, using the workflow config and inputs config to override 
        
        
#            if self.config['use_snakemake_api']:
#                snakemake(self.snakefile,configfiles=[self.updated_config], workdir=config['output_dir'], dryrun=True, debug_dag=True)
#            else:
            #run snakemake command-line (passing any leftover args from argparse)
        snakemake_cmd_list = ['snakemake',
                                f'--snakefile {self.snakefile}',
                                f"--directory {self.config['output_dir']}",
                                f'--configfile {self.updated_config} ',
                                *self.config['snakemake_args'],
                                *self.config['targets_by_analysis_level'][analysis_level]]

        snakemake_cmd = ' '.join(snakemake_cmd_list)
        run(snakemake_cmd)

   
