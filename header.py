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
import itertools 


import scienceplots
import matplotlib.pyplot as plt 

# plt.style.use('science')
# plt.style.use('nature')
plt.style.use('science')
plt.rcParams['image.cmap'] = 'gray'
plt.rcParams["figure.figsize"] =(3.5  , 2.625 )
plt.rcParams['lines.markersize'] = plt.rcParams['lines.markersize'] *0.8
plt.rcParams['lines.linewidth'] /= 1.5



from cycler import cycler
colors=["#ef476f","#f78c6b","#ffd166","#06d6a0","#118ab2","#073b4c"]
ms = [
    '.',   # point
    '<',   # triangle_left
    ',',   # pixel
    'v',   # triangle_down
    'p',   # pentagon
    '^',   # triangle_up
    '>',   # triangle_right
    's',   # square
    'o',   # circle
    '8',   # octagon
]
marker = itertools.cycle(ms)
ls = [':','--','-.','-','.']

cy=cycler('color',colors)*cycler('linestyle', ['-']) 
plt.rcParams['axes.prop_cycle'] = cy  


##functions for mean neaerst distance
from mpmath import gammainc,gamma
def rm(phi):
    return 4**(1/3) * np.sqrt(gammainc(5/3,8*phi))*np.exp(4*phi)*2**(-5/3)*phi**(-1/3)
