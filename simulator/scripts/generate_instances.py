DESC='''Pipeline to generate simulation instances given a config file (See
../input/experiments_test.json) with model parameters and their values.
This is a full factorial design. A folder tree will be created inside a
folder called "./experiments/". Each leaf folder of this tree will contain
the necessary files to run an experiment.'''

import argparse
from configparser import ConfigParser
from glob import glob
from itertools import product
from json import load,dump
import logging
from os import environ
from os.path import basename,exists
from pathlib import Path
from pdb import set_trace
from re import sub

SIMULATOR=None
EXPERIMENTS=None
SUITABILITY_THRESH = 0

FORMAT='[%(asctime)-15s] [%(filename)s] [%(levelname)s]: %(message)s'
PARAMETERS_AND_ATTRIBUTES=[\
        ('number_of_simulations',int), \
        ('exposure_delay',float), \
        ('alpha_S',float), \
        ('alpha_L',float), \
        ('alpha_LD',float), \
        ('moore_range',int), \
        ('start_month',int), \
        ('random_seed',int),\
        ('seeding',str),\
        ]

PARAMETERS=[x[0] for x in PARAMETERS_AND_ATTRIBUTES]
PARAMETER_TYPES=[x[1] for x in PARAMETERS_AND_ATTRIBUTES]
       
LOGFILE='log'
CONFIGFILE='config.json'
OUT_PREFIX = 'out'

def generate_instance(networkPath,instance,slurmFile,dryRun,shell):

    parDict=dict(zip(PARAMETERS,instance))

    #Intialising sub dictionaries
    configDict = {}
    configDict['network_specific_input']={
            "network": networkPath,
            "seeding": parDict['seeding'],
            }
    
    configDict['simulation_parameters'] = {
            "time_steps" : parameters['time_steps'],
            "number_of_simulations" : parDict['number_of_simulations'],
            "start_month" : parDict['start_month'],
            }

    configDict['random_seed'] = parDict['random_seed']
    configDict['model_parameters'] = {
            "suitability_thresh": SUITABILITY_THRESH,
            "exposure_delay": parDict['exposure_delay'],
            "alpha_S": parDict['alpha_S'],
            "alpha_L": parDict['alpha_L'],
            "alpha_LD": parDict['alpha_LD'],
            "range_type": "moore",
            "moore_range": parDict['moore_range']
            }
    configDict['simulation_output_prefix']=OUT_PREFIX

    # create folder
    networkName=sub(parameters['unwanted_network_prefix'],'',networkPath)
    folder=EXPERIMENTS+'/'+'network-'+sub('/','-',networkName)

    for p,v in parDict.items():
        try:
            if p in {'seeding'}:
                v=basename(v)
            if '/' in v:
                v=sub('/','-',v)
            v = sub('.csv','',v)
        except:
            pass
        folder+='/'+p+'-'+str(v) 

    folder=sub(' ','-',folder)

    # create slurm command.
    logFileFullPath=folder+'/'+LOGFILE

    logging.info("Creating %s ..." %folder)
    Path(folder).mkdir(parents=True,exist_ok=True)
    fn = folder+"/"+CONFIGFILE
    
    # Adding the config file to the leaf nodes.
    with open(fn,'w') as fp:
        dump(configDict,fp,indent=4)

    # If outfile exists, skip
    if exists(logFileFullPath):
        with open(logFileFullPath,'r') as lf:
            if 'Done.' in lf.read():
                logging.info("Skipping %s ..." %folder)
                return False
    
    command=f"python {simulatorFolder}/run_spread.py -s config.json"

    if shell:   # invoke as a simple bash command
        slurmFile.write('cd %s; %s > %s\n\n; cd -' %(folder,command,LOGFILE))
    else:
        slurmFile.write(f'''sbatch -D {folder} \
-o {LOGFILE} --export=ALL,command="{command}" \
{simulatorFolder}/run_proc.sbatch; \
{simulatorFolder}/qreg\n''')

    return True

if __name__== "__main__":
    # parser
    parser=argparse.ArgumentParser(description=DESC, 
            formatter_class=argparse.RawTextHelpFormatter)

    # feature vector arguments
    parser.add_argument("config_file",help="Config file")
    parser.add_argument("-r", "--run_file", 
            default="run.sh",
            help="It contains slurm invocation of training script for each instance.")

    parser.add_argument("-n","--dry_run", action="store_true")
    parser.add_argument("--shell", action="store_true",help="Plain shell command, no sbatch.")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args()

    # set logger
    if args.debug:
       logging.basicConfig(level=logging.DEBUG,format=FORMAT)
    elif args.quiet:
       logging.basicConfig(level=logging.WARNING,format=FORMAT)
    else:
       logging.basicConfig(level=logging.INFO,format=FORMAT)

    # check if necessary environment variables are set
    SIMULATOR=environ.get('SIMULATOR')
    EXPERIMENTS=environ.get('EXPERIMENTS')

    if SIMULATOR is None:
        raise ValueError('The environment variable SIMULATOR is not set.')
    elif EXPERIMENTS is None:
        raise ValueError('The environment variable EXPERIMENTS is not set.')

    simulatorFolder=f"{SIMULATOR}/simulator/scripts"

    # read config file in JSON format
    logging.info("Reading JSON file '%s' ..." %args.config_file)
    with open(args.config_file) as f:
        parameters = load(f)

    # create tuples
    ## if parameters['seeding_specification']=='folder':
    ##     parameters['experiments']['seeding']= \
    ##             glob(f"{dataFolder}/{parameters['experiments']['seeding']}/*")
    ## else:
    parameters['experiments']['seeding']= \
            [f"{parameters['experiments']['seeding']}"]

    instances=product(*[parameters['experiments'][k] for k in PARAMETERS])

    # generate run instance for each tuple
    numInstancesProcessed=0
    with open(args.run_file,'w') as f:
        f.write('#!/bin/bash\n')
        f.write('start=$SECONDS\n')
        for instance in instances:
            numInstancesProcessed+=generate_instance(
                    f"{parameters['network']}",
                    instance,f,args.dry_run,args.shell)
        f.write('echo "Total time" $(($SECONDS-$start))\n')
    
    logging.info(f"Number of instances processed: {numInstancesProcessed}")
