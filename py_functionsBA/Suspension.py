from datetime import time
import math
from re import U
from py_functionsBA.header import loadbar
from py_functionsBA.header import timedeco
import matplotlib.pyplot as plt 
import numpy as np
import numpy.linalg as LA
import pandas as pd
from py_functionsBA.Bubble import * 
import os
import sys
import filecmp
class Suspension:
    @timedeco
    def __init__(self,dircase: str,tmin =  0.1,resultsdir ="../results/"):
        print('Reading : '+dircase+'...')
        # check if the studies has already been carried 
        self.resultsdir = resultsdir
        self.PATH = self.resultsdir + dircase +'/'
        # reading of main parameters 
        sys.path.append('./'+self.PATH)
        print(os.listdir(self.PATH))
        import parameters as pstud
        from importlib import reload
        reload(pstud)
        self.parameters = pstud.parameters
        sys.path.remove('./'+self.PATH)
        self.Nb = self.parameters['Nb']      #//int(self.Para['Nb'].values)
        self.L0 = self.parameters['Ls']      #//float(self.Para['Ls'].values)
        self.dt = self.parameters['dtprint']      #//float(self.Para['dt'].values)
        self.g = 1.
        # setting up the data 
        self.tmin = round(tmin,2)
        self.bubbles=list()
        print('reading Pos...')
        self.Pos = pd.read_csv(self.PATH+'Pos.csv');

        # set the mean diameters of bubbles 
        self.D = self.parameters['D']
        self.R = self.parameters['D']/2
        self.dp = self.L0/math.sqrt(self.Nb)
        self.Gs_ad = np.linspace(-1,2,50) # adimensionalised range
        #rearranging Pos tab

        # Reading the Vels file and set all steps
        print('reading Stats...')
        self.Vels = pd.read_csv(self.PATH+'Stats.csv')
        # Reading Dist file
        file = self.PATH+'Dist.csv'
        print('reading Dist...')
        self.Dist = pd.read_csv(file)
        self.Data = self.Dist[self.Dist.t.ge(self.tmin)]
        self.DataRaw = self.Data.drop(columns=['i','t'])
        
        # samples parameters 
        self.TMAX = self.Vels.t.values[-1]
        self.range = 10.
        self.pair_number = self.Nb*(self.Nb-1)/2
        self.PDF,self.PDFx,self.CDF = {},{},{}
    
    def calcul_all(self):
        self.calcul_radial_distrib(Dr=self.D/10)
        self.calcul_contact_frequency()
        self.calcul_contact_time()
        self.calcul_stats_Vels() 
        # self.UprimeUprime_f()
        # self.UprimeUprime_p()
        # self.OmegapOmegap_p()
    
        # self.calcul_relative_vel()
        # self.calcul_acc()
        self.create_bubbles()
        

    def create_bubbles(self):
        for i in range(self.Nb):
            self.bubbles.append(Bubble(i))
            self.bubbles[i].Pos = self.Pos[self.Pos.tag == i].sort_values(['t']).copy()
        
    def calcul_acc(self):
        self.Pos = self.Pos.sort_values(['tag','t'])
        self.Pos['ax'] = self.Pos['vx'].diff()/self.dt
        self.Pos['ax_rel'] = self.Pos['vx_rel'].diff()/self.dt
        self.Pos['ay'] = self.Pos['vy'].diff()/self.dt
        self.Pos['ay_rel'] = self.Pos['vy_rel'].diff()/self.dt
        self.Pos['aomegaz'] = self.Pos['omegaz'].diff()/self.dt
        self.Pos = self.Pos.sort_values(['t','tag'])
        self.Pos = self.Pos[self.Pos.t != 0.1]
        
    def calcul_relative_vel(self):
        self.Pos = self.Pos.sort_values(['t','tag'])
        dtinv = int(1./self.dt)
        self.Pos['tagminall'] = ( self.Pos.tagmin + (self.Pos.t - self.dt) * dtinv * self.Nb)
        self.Pos = self.Pos.astype({"tagminall": int})
        self.Pos['x_rel'] = self.Pos.x.values - self.Pos.iloc[self.Pos.tagminall].x.values
        self.Pos['y_rel'] = self.Pos.y.values - self.Pos.iloc[self.Pos.tagminall].y.values
        self.Pos['vx_rel'] = self.Pos.vx.values - self.Pos.iloc[self.Pos.tagminall].vx.values
        self.Pos['vy_rel'] = self.Pos.vy.values - self.Pos.iloc[self.Pos.tagminall].vy.values
        self.Pos['omegaz_rel'] = self.Pos.omegaz.values - self.Pos.iloc[self.Pos.tagminall].vy.values
             
    @timedeco
    def collision_corr(self):
        # data[[index for index in data.keys() if '1' in index.split('b')]]
        return 1


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
            self.Hz.append(Nc /Dt  * math.sqrt(self.D/self.parameters['g']))
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
        self.Stats = self.Vels[self.Vels.t.ge(self.tmin)].drop(columns=['i','t']).describe()
        self.Stats.to_csv(self.PATH+'meanstats.csv')
        # valocities fluctuation 
     
        self.Vels['dVy'] = self.Vels['vydrops']-self.Vels['vyfluid']
        self.Vels['dVx'] = self.Vels['vxdrops']-self.Vels['vxfluid']
    
    def calcul_PDF(self,list,range = None):
        bins = 100
        H,X1=np.histogram(list,bins=bins,range=range)
        dx = X1[1]-X1[0]
        CDF = np.cumsum(H)/len(list)
        PDFx,PDF = X1[1:],np.gradient(CDF,dx)
        return PDFx,CDF,PDF

    def UprimeUprime_p(self):
        """calculation of the particular average pseudo turbulent tensor
        """
        #time averaged mean

        self.Pos['Uprime_x'] = self.Pos['vx'] - self.Pos['vx'].mean()
        self.Pos['Uprime_y'] = self.Pos['vy'] - self.Pos['vy'].mean()
        # Uprime_p = zip(self.Pos['Uprime_x'].values,self.Pos['Uprime_y'].values)
        # outer = np.array((2,2))
        # for vec in Uprime_p:
        #     outer = np.add(outer,np.outer(vec,vec))
        self.Pos['UU_pxx'] = self.Pos['Uprime_x']**2
        self.Pos['UU_pyy'] = self.Pos['Uprime_y']**2
        self.Pos['UU_pxy'] = self.Pos['Uprime_x']*self.Pos['Uprime_y']
        
        # PDF of the fluct
        self.PDFx['Uprimex'],self.CDF['Uprimex'],self.PDF['Uprimex'] = self.calcul_PDF(self.Pos['Uprime_x'].values)
        self.PDFx['Uprimey'],self.CDF['Uprimey'],self.PDF['Uprimey'] = self.calcul_PDF(self.Pos['Uprime_y'].values)

    def OmegapOmegap_p(self):
        """calculation of the particular average pseudo turbulent tensor
        """
        #time averaged mean

        self.Pos['Omegap'] = self.Pos['omegaz'] - self.Pos['omegaz'].mean()
        self.Pos['OO_pzz'] = self.Pos['Omegap']**2
        
        # PDF of the fluct
        self.PDFx['omegap'],self.CDF['omegap'],self.PDF['omegap'] = self.calcul_PDF(self.Pos['Omegap'].values)

