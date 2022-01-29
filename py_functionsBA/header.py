from typing import Optional
import numpy as np
import pandas as pd
import math
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import filecmp
import time

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



#########################################################
#########    Plots of the suspensions class    ##########
#########################################################
def plot_all_fig(allStudys: list,Alllabel: list, PATH =''):
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
        # plt.hist(suspension_ga10_mu004.time_of_contacts,bins= 50,density=True,label = r'$N = 10$')
    plt.xlabel(r'$t/\sqrt{D/g}$')
    plt.ylabel(r'$PDF$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/timePDF.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.PDFx['flucx'],std.PDF['flucx'],label = label)
    plt.xlabel(r'$U_x$')
    plt.ylabel(r'$PDF$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/velfluctuationX.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.PDFx['flucy'],std.PDF['flucy'],label = label)
    plt.xlabel(r'$U_y$')
    plt.ylabel(r'$PDF$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/velfluctuationY.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['dVy'],label = label)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$<U_y>$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/RelativeVelY.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['dissAir']/std.L0**2,label = label)
    # plt.savefig(PATH+'/velfluctuationtime.pdf', format='pdf',bbox_inches='tight')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$||\overline{\tau_{fluid}}||/L^2$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/taufluid.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['dissWater']/std.L0**2,label = label)
    # plt.savefig(PATH+'/velfluctuationtime.pdf', format='pdf',bbox_inches='tight')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$||\overline{\tau_{bubbles}}||/L^2$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/taububbles.pdf', format='pdf',bbox_inches='tight')
    plt.show()

    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['rhovx'],label = label)
    # plt.savefig(PATH+'/velfluctuationtime.pdf', format='pdf',bbox_inches='tight')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\rho u_x$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/rhovx.pdf', format='pdf',bbox_inches='tight')
    plt.show()

    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.Vels['t'],std.Vels['rhovy'],label = label)
    # plt.savefig(PATH+'/velfluctuationtime.pdf', format='pdf',bbox_inches='tight')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\rho u_y$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/rhovy.pdf', format='pdf',bbox_inches='tight')
    plt.show()
    for std,label in zip(allStudys,Alllabel):
        plt.plot(std.volb,label = label)
    # plt.savefig(PATH+'/velfluctuationtime.pdf', format='pdf',bbox_inches='tight')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$V/V_{ini}$')
    plt.legend()
    if PATH:
        plt.savefig(PATH+'/vol.pdf', format='pdf',bbox_inches='tight')
    plt.show()