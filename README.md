# Code for "Network Models and Simulation Analytics for Multi-scale Dynamics of Biological Invasions, Adiga et al. 2021, Frontiers BigData"
Contact: Abhijin Adiga (abhijin@virginia.edu, abhijin@gmail.com)

The experiments were performed using a high-performance computing cluster in Linux environment (https://www.rc.virginia.edu/userinfo/rivanna/overview/). Therefore, the user will encounter some SLURM commands. These can be easily replaced by normal Unix shell statements.

## Synthetic network generation
Folder: ``./synthetic_networks``
See README.md in this folder for information.

## Multi-pathway simulations
Folder: ``./simulator``
See README.md in this folder for information.

## Structural and dynamical analysis
Folder: ``./analysis``
See README.md in this folder for information.

## Results
The results of network measurements and simulations are in ``results.db``. It has two tables.

* ``mpsn_props`` has network measures for all the generated synthetic networks.
* ``summary`` has simulation results.
