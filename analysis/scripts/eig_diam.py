DESC='''
Plot eigenvalues and diameter
'''

from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from os import environ,getenv
from os.path import basename
import pandas as pd
from pdb import set_trace
import seaborn as sns

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
FORMAT='[%(asctime)-15s] [%(filename)s] [%(levelname)s]: %(message)s'

# read data
df_list = []
for f in glob('../results/eig_diam/*csv'):
    df_list.append(pd.read_csv(f))
mdf = pd.concat(df_list).sort_values('network')

for r in range(1,3):
    df = mdf[mdf.range == r]
    df.network = df.network.apply(basename).str.replace('.csv','')
    df = df.rename(columns = {
        'spectral_radius': 'Weighted spectral radius',
        'spectral_radius_unweighted': 'Unweighted spectral radius',
        'diameter': 'Diameter',
        'network': 'Network',
        'month': 'Month'
        })
    df = df[(df.Network != 'SN') & (df.Network != 'NP')]
    
    # weighted spectral radius
    sns.lineplot(data = df, x = 'Month', y = 'Weighted spectral radius', 
            hue = 'Network', markers = True)
    plt.title(f'Range {r}')
    plt.savefig(f'real_weighted_spectral_radius_{r}.pdf', bbox_inches = 'tight')
    plt.clf()
    
    # weighted spectral radius
    sns.lineplot(data = df, x = 'Month', y = 'Unweighted spectral radius', 
            hue = 'Network', markers = True)
    plt.title(f'Range {r}')
    plt.savefig(f'real_unweighted_spectral_radius_{r}.pdf', bbox_inches = 'tight')
    plt.clf()
    
    # diameter
    df = df[df.Diameter != -1]
    sns.lineplot(data = df, x = 'Month', y = 'Diameter', 
            hue = 'Network', markers = True)
    plt.title(f'Range {r}')
    plt.savefig(f'real_diameter_{r}.pdf', bbox_inches = 'tight')
    plt.clf()
