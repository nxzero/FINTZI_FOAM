import pandas as pd
import numpy as np
import sympy as sp
from numpy import linalg as LA
from scipy import optimize
from scipy.optimize import curve_fit
import itertools
import math
import seaborn as sns
from matplotlib.colors import PowerNorm
import os
import sys
from mpmath import gammainc,gamma

# mp.rcParams.update({
#     'text.usetex' : True,
#     'font.family': 'sans-serif',
#     'font.sans-serif': ['Helvetica'],
#     'font.size': 10,
#     'axes.grid': True,
#     'grid.linestyle': '--',
#     'grid.linewidth': 0.5,
#     'grid.alpha': 0.5,
#     'xtick.major.width': 0.5,
#     'xtick.minor.width': 0.25,
#     'ytick.major.width': 0.5,
#     'ytick.minor.width': 0.25,
#     'legend.frameon': False,
#     'figure.figsize': [3.39, 2.1],
#     'lines.markersize': 5,
#     'lines.linewidth':1,
# })
# mp.rcParams['savefig.format'] = 'pdf'
# mp.rcParams['savefig.bbox'] = 'tight'
import scienceplots
import matplotlib.pyplot as plt 

# plt.style.use('science')
# plt.style.use('nature')
plt.style.use('science')
plt.rcParams['image.cmap'] = 'turbo'
plt.rcParams["figure.figsize"] =(3.5  * 0.8, 2.625  * .8)
plt.rcParams['lines.markersize'] = plt.rcParams['lines.markersize'] *0.8
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
ms = ['o','^','s','d','>','<']
ls = [':','--','-.','-','.']

##functions for mean neaerst distance
from mpmath import gammainc,gamma
def rm(phi):
    return 4**(1/3) * np.sqrt(gammainc(5/3,8*phi))*np.exp(4*phi)*2**(-5/3)*phi**(-1/3)
