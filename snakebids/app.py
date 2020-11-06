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





class SnakeBidsApp:

    def __init__(self, snakebids_config=None, snakefile=None):
    
        if snakebids_config == None:
            #try to locate config file relative to this file 
            this_dir = os.path.dirname(os.path.realpath(__file__))
            print(f'this dir is {this_dir}')
            self.snakebids_config = os.path.join(this_dir,'config','snakebids.yml')
        else:
            self.snakebids_config = snakebids_config

        #make sure it exists
        if not os.path.exists(self.snakebids_config):
            raise FileNotFoundError('snakebids_config file not found must be defined or placed in config/snakebids.yml')


        if snakefile == None:
            #try to locate config file relative to this file
            this_dir = os.path.dirname(os.path.realpath(__file__))

            #could replace this with code that looks in specific locations (e.g. from snakemake)
            self.snakefile = os.path.join(this_dir,'workflow','Snakefile')
        else:
            self.snakefile = snakefile

        #make sure snakefile exists
        if not os.path.exists(self.snakefile):
            raise FileNotFoundError('Snakefile not found must be defined or placed in workflow/Snakefile')

        self.parse_args()
        self.write_updated_config()


    def parse_args(self):

        #replace with proper name of pipeline here
        parser = argparse.ArgumentParser(description='snakebids-app')

        #load up workflow config file
        with open(self.snakebids_config, 'r') as infile:
            self.config = yaml.load(infile, Loader=yaml.FullLoader)

        #update the parser with config options
        for name, parse_args in self.config['parse_args'].items():
            parser.add_argument(name, **parse_args)

        all_args = parser.parse_known_args()

        args = all_args[0]
        snakemake_args = all_args[1]

        #add arguments to config
        self.config.update(args.__dict__)
        self.config.update({'snakemake_args': snakemake_args})

        

    def write_updated_config(self):
        #create an updated snakebids config file
        self.updated_config = os.path.join(self.config['output_dir'],'config','snakebids.yml')

        #write it to file
        os.makedirs(os.path.dirname(self.updated_config),exist_ok=True)
        with open(self.updated_config, 'w') as outfile:
            yaml.dump(self.config, outfile, default_flow_style=False)




    # run snakemake with that config
    #workflow snakefile will read snakebids config, create inputs_config, read that in
    def run_snakemake(self):
        
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

   
