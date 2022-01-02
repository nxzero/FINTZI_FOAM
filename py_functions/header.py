import os
import shutil
import re
import math
from PP_class import *
from moment_Cox import *
from triPP_class import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import mpmath
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

import numpy  as np 
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scipy.special as sc
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.size":30,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Helvetica"]})
# # for Palatino and other serif fonts use:
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.size": 30,
#     # "font.family": "serif",
#     "font.serif": ["Helvetica"],
#     "lines.linewidth" : 3,
#     "lines.markersize" : 10,
    # "legend.fontsize":20,
    # 'axes.color_cycle':['red','#fb8500','#219ebc','#023047','#8ecae6','#ffb703'],
    # 'axes.color_cycle':['#0072bd','#d95319','#edb120','#7e2f8e','#77ac30','#4dbeee','#a2142f']
    # 'axes.color_cycle':["#001219","#005f73","#0a9396","#94d2bd","#e9d8a6","#ee9b00","#ca6702","#bb3e03","#ae2012","#9b2226"],
    # 'axes.color_cycle':["#ef476f","#ffd166","#06d6a0","#118ab2","#073b4c"],
    # 'axes.color_cycle':["#d61900","#ff9d2e","#ffd042","#002db3","#1f87ff","#00b0c7","#1cd100","#800094"]
# })
plt.rcParams.update({
    "font.size":30,
    "lines.linewidth" : 3,
    "lines.markersize" : 10,
    "text.usetex": True,
    "font.family": "Helvetica"
})
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["#d61900","#ff9d2e","#ffd042","#002db3","#1f87ff","#00b0c7","#1cd100","#800094"]) 



Thetas = [0,10,20,30,40,50,60,70,80,90]
Chis = [2,5,10,15,20,30,40]
S_start2 ={}
Ress = []
Thetass = []
for Chi in Chis:
    if Chi == 40:
        Thetas = [40,50]
    Thetass.append(Thetas)
    S_start2[str(Chi)] = {}
    for Theta in Thetas:
        S_start2[str(Chi)][str(Theta)] = Study('RTC/ReVsThetha_Chi'+str(Chi)+'/Re_Theta'+str(Theta)+'/')
        S_start2[str(Chi)][str(Theta)].get_forces_and_troques()
    if 0.005 in S_start2[str(Chi)][str(Theta)].parameters_list:
        Ress.append(S_start2[str(Chi)][str(Theta)].parameters_list[1:])#
    else:
        Ress.append(S_start2[str(Chi)][str(Theta)].parameters_list)#
S2 = {}
for Chi,Res,Thetas in zip(Chis,Ress,Thetass):
    S2[str(Chi)] = {}
    for Re in Res:
        Theta = Thetas[0]
        i = S_start2[str(Chi)][str(Theta)].parameters_list.index(Re)
        S_temp_para  = S_start2[str(Chi)][str(Theta)].parameters_list[i]
        S2[str(Chi)][str(S_temp_para)] = {}
        S2[str(Chi)][str(S_temp_para)]['parameters_list'] = []
        S2[str(Chi)][str(S_temp_para)]['F_dpara'] = []
        S2[str(Chi)][str(S_temp_para)]['F_d2para'] = []
        S2[str(Chi)][str(S_temp_para)]['F_dperp'] = []
        S2[str(Chi)][str(S_temp_para)]['F_d2perp'] = []
        S2[str(Chi)][str(S_temp_para)]['M_had'] = []
        S2[str(Chi)][str(S_temp_para)]['eM_had'] = []
        S2[str(Chi)][str(S_temp_para)]['F_d'] = []
        S2[str(Chi)][str(S_temp_para)]['eF_d'] = []
        S2[str(Chi)][str(S_temp_para)]['F_dd'] = []
        S2[str(Chi)][str(S_temp_para)]['F_l'] = []
        S2[str(Chi)][str(S_temp_para)]['eF_l'] = []
        S2[str(Chi)][str(S_temp_para)]['F_dl'] = []
    for Theta in Thetas:
        for Re in Res:
            i = S_start2[str(Chi)][str(Theta)].parameters_list.index(Re)
            S_temp_para  = S_start2[str(Chi)][str(Theta)].parameters_list[i]
            S_temp_F_dpara =S_start2[str(Chi)][str(Theta)].F_dpara[i]
            S_temp_F_d2para =S_start2[str(Chi)][str(Theta)].F_d2para[i]
            S_temp_F_dperp  =S_start2[str(Chi)][str(Theta)].F_dperp[i]
            S_temp_F_d2perp  =S_start2[str(Chi)][str(Theta)].F_d2perp[i]
            S_temp_M_had  =S_start2[str(Chi)][str(Theta)].M_had[i]
            S_temp_eM_had  =S_start2[str(Chi)][str(Theta)].eM_had[i]
            S_temp_F_d  =S_start2[str(Chi)][str(Theta)].F_d[i]
            S_temp_eF_d  =S_start2[str(Chi)][str(Theta)].eF_d[i]
            S_temp_F_dd  =S_start2[str(Chi)][str(Theta)].F_dd[i]
            S_temp_F_l  =S_start2[str(Chi)][str(Theta)].F_l[i]
            S_temp_eF_l  =S_start2[str(Chi)][str(Theta)].eF_l[i]
            S_temp_F_dl  =S_start2[str(Chi)][str(Theta)].F_dl[i]
            S2[str(Chi)][str(S_temp_para)]['parameters_list'].append(Theta)
            S2[str(Chi)][str(S_temp_para)]['F_dpara'].append(S_temp_F_dpara)
            S2[str(Chi)][str(S_temp_para)]['F_d2para'].append(S_temp_F_d2para)
            S2[str(Chi)][str(S_temp_para)]['F_dperp'].append(S_temp_F_dperp)
            S2[str(Chi)][str(S_temp_para)]['F_d2perp'].append(S_temp_F_d2perp)
            S2[str(Chi)][str(S_temp_para)]['M_had'].append(S_temp_M_had)
            S2[str(Chi)][str(S_temp_para)]['eM_had'].append(S_temp_eM_had)
            S2[str(Chi)][str(S_temp_para)]['F_d'].append(S_temp_F_d)
            S2[str(Chi)][str(S_temp_para)]['eF_d'].append(S_temp_eF_d)
            S2[str(Chi)][str(S_temp_para)]['F_dd'].append(S_temp_F_dd)
            S2[str(Chi)][str(S_temp_para)]['F_l'].append(S_temp_F_l)
            S2[str(Chi)][str(S_temp_para)]['eF_l'].append(S_temp_eF_l)
            S2[str(Chi)][str(S_temp_para)]['F_dl'].append(S_temp_F_dl)



markers = [ 'D', 'H','s','d','^','<','o', 'p','h', 'v', '>','*','8', 'd']
# colors = ["#d61900","#00b0c7","#ff7c6b","#ff9d2e","#ffd042","#002db3","#1f87ff","#1cd100","#800094"]
colors=["#d61900","#ff9d2e","#ffd042","#002db3","#1f87ff","#00b0c7","#1cd100","#800094"]

