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
    def __init__(self,dircase: str,tmin =  0.1,resultsdir ="../results/" ):
        print('Reading : '+dircase+'...')
        # check if the studies has already been carried 
        self.resultsdir = resultsdir
        self.PATH = self.resultsdir + dircase +'/'
        file = self.PATH+'Pos.csv'
        parafile = self.PATH+'Para.csv'
        # reading of main parameters 
        self.Para = pd.read_csv(parafile)
        self.L0 = float(self.Para['Ls'].values)
        self.dt = float(self.Para['dt'].values)
        self.g = 1.
        self.Pos = pd.read_csv(file)
        # setting up the data 
        self.TimeList = list(set(self.Pos.t.values))
        self.TimeList.sort()
        self.tmin = round(tmin,2)
        index = self.TimeList.index(tmin)
        self.TimeList = self.TimeList[index:]
        self.ListStep = list()
        self.Pos = self.Pos[self.Pos.t.ge(self.tmin)]
        # separate all the steps 
        for t in self.TimeList:
            self.ListStep.append(self.Pos[self.Pos.t.eq(t)])

        self.tagmax = max(self.Pos.tag.values)
        self.bubbles=list()
        for tag in range(self.tagmax+1):
            b = Bubble(tag)
            bstate = self.ListStep[0][self.ListStep[0].tag.eq(tag)]
            b.assigne(bstate)
            self.bubbles.append(b)
            
        # set the mean diameters of bubbles 
        self.D = float(self.Para['D'].values)
        self.Nb = len(self.bubbles)
        self.dp = self.L0/math.sqrt(self.Nb)
        self.Gs_ad = np.linspace(-1,2,50) # adimensionalised range
        # Reading the Vels file and set all steps
        
        self.pair_list = []
        for i in range(self.tagmax):
            for j in range(i+1,self.tagmax+1):
                self.pair_list.append('b'+str(i)+'b'+str(j))
        self.pair_number = len(self.pair_list)
        self.Vels = pd.read_csv(self.PATH+'Stats.csv')
        # Reading Dist file
        file = self.PATH+'Dist.csv'
        self.Dist = pd.read_csv(file)
        self.Data = self.Dist[self.Dist.t.ge(self.tmin)]
        self.DataRaw = self.Data.drop(columns=['i','t'])
        
        # samples parameters 
        self.TMAX = self.Vels.t.values[-1]
        self.range = 10.
        
        self.PDF,self.PDFx,self.CDF = {},{},{}
    
    def calcul_all(self):
        self.create_list()
        self.calcul_radial_distrib(Dr=self.D/10)
        self.calcul_contact_frequency()
        self.calcul_contact_time()
        self.calcul_stats_Vels() 

    
    @timedeco
    def create_list(self):
        """reads the Pos.csv file and create the list of bubbles with the rigth pos.

        Args:
            plot (bool, optional): [Plot the path of bubbles into a pdf file called trajet.pdf]. Defaults to False.
        """
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
        # self.Vels['dVx'] = self.Vels['dissAir']
        # self.Vels['dVx'] = self.Vels['dissWater']
    
    def calcul_PDF(self,list,range = None):
        bins = 100
        H,X1=np.histogram(list,bins=bins,range=range)
        dx = X1[1]-X1[0]
        CDF = np.cumsum(H)/len(list)
        PDFx,PDF = X1[1:],np.gradient(CDF,dx)
        return PDFx,CDF,PDF

