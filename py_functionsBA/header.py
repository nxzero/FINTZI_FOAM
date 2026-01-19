from typing import Optional
import numpy as np
import pandas as pd
import math
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import filecmp
import time
import matplotlib as mpl
RESULT_DIR='/home/irsrvshare1/R04/RheoPipe/fintzin/CYLINDERS_PROJECT.backup/results/'

def loadbar(iterration, total, prefix='',suffix='',decimals=1, length=50,fill='>'):
    if iterration == 0:
        percent=('{0:.'+str(decimals)+'f}').format(100* (iterration/float(total)))
        filledLength = int(length*iterration//total)
        bar = fill * filledLength +'-'*(length-filledLength)
        print(f'\r{prefix}|{bar}|0%{suffix}',end='\r')
    else:
        percent=('{0:.'+str(decimals)+'f}').format(100* (iterration/float(total)))
        filledLength = int(length*iterration//total)
        bar = fill * filledLength +'-'*(length-filledLength)
        print(f'\r{prefix}|{bar}|{percent}%{suffix}',end='\r')
    if iterration == total-1:
        percent=('{0:.'+str(decimals)+'f}').format(100* (iterration/float(total)))
        filledLength = int(length*iterration//total)
        bar = fill * filledLength +'>'*(length-filledLength)
        print(f'\r{prefix}|{bar}|100%{suffix}',end='\r')
        print('\n')
        
def timedeco(func):
    def inner(*args, **kwargs):
        default = {'Time':True}
        default.update(kwargs)
        if default['Time']:start = time.time()
        func(*args, **kwargs)
        if default['Time']:stop = time.time();print('Duration :',stop-start,'s')
    return inner

plt.rcParams.update({
    "font.size":18,
    "lines.linewidth" : 2,
    "lines.markersize" : 10,
    "text.usetex": True,
    "font.family": "serif"
})
colors=["#d61900","#ff9d2e","#ffd042","#1cd100","#00b0c7","#1f87ff","#002db3","#800094"]#["#d61900","#ff9d2e","#ffd042","#002db3","#1f87ff","#00b0c7","#1cd100","#800094"]
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors) 


#########################################################
#########    Plots of the suspensions class    ##########
#########################################################
def plot_all_fig(allStudys: list,Alllabel: list, PATH ='',bins = 30):
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Pr,std.P,label = label)
    plt.xlabel(r'$r/D$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/dist.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Gs_ad,std.Hz,label = label)
    plt.xlabel(r'$\frac{G-D}{dp-D}$')
    plt.ylabel(r'$\frac{f}{V \Delta t}$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/Hz.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.PDFx['ct'],std.PDF['ct'],label = label)
    plt.xlabel(r'$t/\sqrt{D/g}$')
    plt.ylabel(r'$PDF$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/timePDF.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
    #Velocities fluctuations 
    # Up_XX**2
    
    
    
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['dissF']/std.L0**2,label = label)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$||\overline{\tau_{fluid}}||/L^2$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/taufluid.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['dissD']/std.L0**2,label = label)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$||\overline{\tau_{bubbles}}||/L^2$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/taububbles.pdf', format='pdf',bbox_inches='tight')
    plt.show()

    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['rhovx'],label = label)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\rho u_x$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/rhovx.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['rhovy'],label = label)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\rho u_y$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/rhovy.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    

    # for std,label in zip(allStudys,Alllabel):
    #     plt.plot(std.volb,label = label)
    # plt.xlabel(r'$t$')
    # plt.ylabel(r'$V/V_{ini}$')
    # plt.legend()
    # if PATH:
    #     plt.savefig(PATH+'/vol.pdf', format='pdf',bbox_inches='tight')
    # plt.show()
    for std,label in zip(allStudys,Alllabel):
        plt.hist(std.Vels['t'],std.Vels['Vx'],label = label)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$<u>$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/Vx.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['Vy'],label = label)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$<u>$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/Vy.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
def color_negative_red(val):
    if abs(val) < 0.2:  
        color = 'black' 
    elif abs(val) < 0.6: 
        color = 'blue' 
    elif abs(val) < 0.9: 
        color = 'orange' 
    elif abs(val) > 0.9: 
        color = 'red' 
    else:
        color = 'green' 
        
    return 'color: %s' % color


# t = s.bubbles[0].Pos['t'].values[750:]
# y = s.bubbles[0].Pos['defnorm'].values[750:]
def FFT(x,y):
    sp = np.abs(np.fft.fft(y))
    freq = np.fft.fftfreq(x.shape[-1],1/(x[0]-x[1])) 
    i = freq > 0
    return freq[i],sp[i]
# plt.plot(freq[i], sp[i],'.')
# plt.xlim([0,0.01])
# plt.show()
