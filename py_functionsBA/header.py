from typing import Optional
import numpy as np
import pandas as pd
import math
import scipy as sp
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
        default = {'Time':False}
        default.update(kwargs)
        if default['Time']:start = time.time()
        func(*args, **kwargs)
        if default['Time']:stop = time.time();print('\nDuration :',stop-start,'s')
    return inner