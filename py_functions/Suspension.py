from datetime import time
import math
from py_functionsBA.header import loadbar
from py_functionsBA.header import timedeco
import matplotlib.pyplot as plt 
import numpy as np
import numpy.linalg as LA
import pandas as pd
from py_functionsBA.Bubble import * 
import os
import filecmp
class Suspension:
    @timedeco
    def __init__(self,dircase: str,tmin =  0.1,resultsdir ="../results/"  ):
        print('Reading : '+dircase+'...')
        # check if the studies has already been carried 
        self.resultsdir = resultsdir
        self.PATH = self.resultsdir + dircase +'/'
        
    @timedeco
    def create_list(self):
        for i,b in enumerate(self.bubbles):
            b.assigne_list(self.Pos)
            loadbar(i, len(self.bubbles), prefix='Pos')
            
            


    @timedeco
    def calcul_radial_distrib(self,Dr = 0.01):
        rmin = self.D/100 # start from 1 pourcent of the radius
        rmax = 3#self.L0/8
        self.P = [] # probability density function 
        self.r = np.arange(rmin,rmax,Dr)
        self.Pr = np.arange(rmin,rmax,Dr)[:-1]/self.D
        number_of_samples = len(self.DataRaw.index) 
        density_of_particles = 1/self.L0**2
        for i,r,rdr in zip(range(len(self.Pr)),self.r[:-1],self.r[1:]):
            sum = ((self.DataRaw>r) & (self.DataRaw<rdr)).sum().sum()
            area = math.pi*(rdr**2-r**2)
            self.P.append(sum/(self.pair_number*area*number_of_samples*density_of_particles))
            loadbar(i+1, len(self.r), prefix='Dist')
    
    @timedeco
    def calcul_contact_frequency(self):
        self.Gs = self.Gs_ad * (self.dp-self.D)+self.D
        self.Hz = []
        for i,G in enumerate(self.Gs):
            DataB = self.DataRaw < self.D * G
            Nc=((DataB != DataB.shift(1)) & DataB).sum().sum()# number of contact
            Dt = self.Data.t.values[-1] - self.tmin
            vol = self.L0**2
            self.Hz.append(Nc /Dt/vol)
            loadbar(i+1, len(self.Gs), prefix='Freq')
    
    def calcul_contact_time(self,G=1):
        DataB = self.DataRaw < self.D * G
        DataP = (DataB != DataB.shift(1))[1:]
        self.DataP =DataP
        self.DataB =DataB
        allsteps = []
        for j,i in enumerate(DataP):
            steps = list(np.where(DataP[i])[0])
            if DataB[i].values[0] == True : steps = steps[1:] # cut first colli 
            if DataB[i].values[-1] == True : steps = steps[:-1] # cut last colli 
            if len(steps) % 2 != 0:
                print('\n !!! !Unpair number of contact\n')
            loadbar(j+1, len(DataP.columns), prefix='CT')
            
            allsteps += steps
        step_start = [x for i,x in enumerate(allsteps) if i % 2 == 0] #gather all pair indicies 
        step_stop = [x for i,x in enumerate(allsteps) if i % 2 == 1] #gather all unpair indicies
        self.time_of_contacts = (np.array(step_stop)-np.array(step_start))*self.dt/math.sqrt(self.D/self.g)
        #self.time_of_contacts=self.time_of_contacts[self.time_of_contacts<100]#in case of bugs maybe
        # density list
        self.PDFx['ct'],self.CDF['ct'],self.PDF['ct'] = self.calcul_PDF(self.time_of_contacts)

    def calcul_stats_Vels(self):
        self.Stats = self.Vels[self.Vels.t.ge(self.tmin)].drop(columns=['i','t','N_Vof']).describe()
        self.Stats.to_csv(self.PATH+'meanstats.csv')
        # valocities fluctuation 
        self.U_fluctx  = self.Pos.vx.values - self.Pos.vx.values.mean()
        self.U_flucty  = self.Pos.vy.values - self.Pos.vy.values.mean()
        # PDF of the fluct
        self.PDFx['flucx'],self.CDF['flucx'],self.PDF['flucx'] = self.calcul_PDF(self.U_fluctx)
        self.PDFx['flucy'],self.CDF['flucy'],self.PDF['flucy'] = self.calcul_PDF(self.U_flucty)
        self.Vels['dVy'] = self.Vels['vydrops']-self.Vels['vyfluid']
        self.Vels['dVx'] = self.Vels['vxdrops']-self.Vels['vxfluid']
        vol = np.zeros(len(self.bubbles[0].vol))
        for b in self.bubbles:
            vol = vol + b.vol
        self.volb = vol/vol[0]
    
    def calcul_PDF(self,list,range = None):
        bins = 100
        H,X1=np.histogram(list,bins=bins,range=range)
        dx = X1[1]-X1[0]
        CDF = np.cumsum(H)/len(list)
        PDFx,PDF = X1[1:],np.gradient(CDF,dx)
        return PDFx,CDF,PDF

