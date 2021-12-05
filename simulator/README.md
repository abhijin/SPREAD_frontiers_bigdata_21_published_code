# Multi-pathway simulator
All code is in ``./scripts``.
``run_spread.py`` is the simulator script which requires the script ``msc_network.py`` to read the networks.

## Setting the environment variables 
``../setenv/frontiers.sh`` can be used as a guide to set the environment variables.

## Generate simulation instances
* Experiments are specified using a config file. See ``input/test_design.json``.
* ``generate_instances.py`` is used to generate the necessary script to run experiments. It currently supports format suitable for the cluster (Slurm commands). However, it can be easily modified for a non-cluster output.
* ``master_experiments.sh test`` has a test instance to illustrate the usage of this framework. 
* The output of ``generate_instances.py`` is 
    * a file called ``run.sh``, a bash script containing a list of commands, each corresponding one simulation configuration. 
    * a folder called ``experiments/`` is created with a nested folders structure. Each leaf folder in this structure corresponds to one simulation configuration specified in the configuration file. 
    * ``bash run.sh`` will run each experiment specified in the leaf folders.
* Outputs can be collected using ``master_experiments.sh collect_results``. The output of this command is ``to_db.sqlite``. Each line is a database statement. The contents of this file can be imported into a database like ``../results.db`` for further analysis..
