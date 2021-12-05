DESC='''
Given scenario, plot number of infections w.r.t. combinations of probabilities.
'''

import argparse
import logging
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from os import environ,getenv
from os.path import basename
import pandas as pd
from pdb import set_trace
import seaborn as sns
import sqlite3

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
FORMAT='[%(asctime)-15s] [%(filename)s] [%(levelname)s]: %(message)s'
COLOR = sns.color_palette()
MARKER = ["o", "^", "s", "v"]
MATHEMATICA=['#5e82b5','#e09c24','#8fb030','#eb634f','#8778b3','#c46e1a','#5c9ec7','#fdbf6f']

def alpha_s(df, prob_list):
    initial_inf = df[(df.alpha_S == 0) &
            (df.alpha_L == 0) &
            (df.alpha_LD == 0)].infections_mean.values[0]
    # base
    p0 = prob_list[1]
    for ind, p in enumerate(prob_list):

        to_plot = df[(df.alpha_S == p) & (df.alpha_L == p0) & (df.alpha_LD == p0)]
        to_plot = to_plot.append({'infections_mean': initial_inf,
            'time_step': 0}, ignore_index=True)
        ## x = [0] + to_plot.time_step.tolist()
        ## y = [initial_inf] + to_plot.infections_mean.tolist()
        ax = sns.lineplot(data = to_plot,
                x = 'time_step', 
                y = 'infections_mean', 
                label = "$\\alpha_S=%g, \\alpha_L=%g, \\alpha_{LD}=%g$" %(p,p0,p0),
                color = MATHEMATICA[0],
                marker = MARKER[ind])

def alpha_l(df, prob_list):
    initial_inf = df[(df.alpha_S == 0) &
            (df.alpha_L == 0) &
            (df.alpha_LD == 0)].infections_mean.values[0]
    # base
    p0 = prob_list[1]
    for ind, p in enumerate(prob_list):

        to_plot = df[(df.alpha_S == p0) & (df.alpha_L == p) & (df.alpha_LD == p0)]
        to_plot = to_plot.append({'infections_mean': initial_inf,
            'time_step': 0}, ignore_index=True)
        ## x = [0] + to_plot.time_step.tolist()
        ## y = [initial_inf] + to_plot.infections_mean.tolist()
        ax = sns.lineplot(data = to_plot,
                x = 'time_step', 
                y = 'infections_mean', 
                label = "$\\alpha_S=%g, \\alpha_L=%g, \\alpha_{LD}=%g$" %(p0,p,p0),
                color = MATHEMATICA[1],
                marker = MARKER[ind])

def alpha_ld(df, prob_list):
    initial_inf = df[(df.alpha_S == 0) &
            (df.alpha_L == 0) &
            (df.alpha_LD == 0)].infections_mean.values[0]
    # base
    p0 = prob_list[1]
    for ind, p in enumerate(prob_list):

        to_plot = df[(df.alpha_S == p0) & (df.alpha_L == p0) & (df.alpha_LD == p)]
        to_plot = to_plot.append({'infections_mean': initial_inf,
            'time_step': 0}, ignore_index=True)
        ## x = [0] + to_plot.time_step.tolist()
        ## y = [initial_inf] + to_plot.infections_mean.tolist()
        ax = sns.lineplot(data = to_plot,
                x = 'time_step', 
                y = 'infections_mean', 
                label = "$\\alpha_S=%g, \\alpha_L=%g, \\alpha_{LD}=%g$" %(p0,p0,p),
                color = MATHEMATICA[2],
                marker = MARKER[ind])

def main():
    # parser
    parser=argparse.ArgumentParser(description=DESC, 
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--locality_size", type = int, required = True, help="Locality size")
    parser.add_argument("--epsilon", type = float, required = True, help="epsilon for probability")
    parser.add_argument("--range", type = int, required = True, help="Range")
    ## parser.add_argument("--time_step", type = int, required = True, help="Time step")
    #parser.add_argument("--seeding", required = True, help="Seeding type: all/loc")
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

    # read table
    conn = sqlite3.connect(getenv('RESULTS_DB'))
    df = pd.read_sql('SELECT * FROM summary',conn)
    df = df[['network',
        'moore_range',
        'alpha_S',
        'alpha_L',
        'alpha_LD',
        'seeding',
        'time_step',
        'infections_mean', 'infections_std',
       'infections_min', 'infections_25_per', 'infections_50_per',
       'infections_75_per', 'infections_max']]
    df['locality_size'] = df.network.str.replace('.*s_','').str.replace('/.*','').astype(int)
    df['epsilon'] = df.network.str.replace('.*p_','').str.replace('/.*','').astype(float)
    ##         (df.time_step == args.time_step)]
    df.seeding = df.seeding.apply(basename).str.replace('.csv','').str.replace('seed_','')
    alpha_s(df[
        (df.locality_size == args.locality_size) & 
        (df.epsilon == args.epsilon) &
        (df.moore_range == args.range)
        ], [0, 0.001, 0.005, 0.01])
    
    alpha_l(df[
        (df.locality_size == args.locality_size) & 
        (df.epsilon == args.epsilon) &
        (df.moore_range == args.range)
        ], [0, 0.001, 0.005, 0.01])

    alpha_ld(df[
        (df.locality_size == args.locality_size) & 
        (df.epsilon == args.epsilon) &
        (df.moore_range == args.range)
        ], [0, 0.001, 0.005, 0.01])

    plt.title(f'locality size $s={args.locality_size}$; $\epsilon= {args.epsilon}$, range $r={args.range}$')
    plt.xlabel('Time')
    plt.ylabel('Average accumulated infections')
    plt.ylim(bottom = 200, top = 500)
    plt.savefig(
            f'non_uniform_s_{args.locality_size}_eps{args.epsilon}_range{args.range}.pdf',
            bbox_inches = 'tight')
    
if __name__== "__main__":
    main()
