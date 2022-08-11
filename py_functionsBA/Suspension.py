import math
from py_functionsBA.header import loadbar
from py_functionsBA.header import timedeco
import numpy as np
# import pandas as pd
import dask.dataframe as pd
import matplotlib.pyplot as plt
# import dask.dataframe as pd
from py_functionsBA.Bubble import * 
from dask import delayed
import sys
import os

class Suspension:
    @timedeco
    def __init__(self,dircase: str,tmin =  0.1,DIST = True,RD =True,resultsdir ="../results/"):
        print('Reading : '+dircase+'...')
        self.DIST = DIST
        self.RD = RD
        # check if the studies has already been carried 
        self.resultsdir = resultsdir
        self.PATH = self.resultsdir + dircase +'/'
        # reading of main parameters 
        sys.path.insert(0,'./'+self.PATH)
        # print(os.listdir(self.PATH))
        import parameters as pstud
        from importlib import reload
        reload(pstud)
        sys.path.remove('./'+self.PATH)
        self.parameters = pstud.parameters
        self.Nb = self.parameters['Nb']      #//int(self.Para['Nb'].values)
        self.rho_f = self.parameters['rho_f']      #//int(self.Para['Nb'].values)
        self.rho_d = self.parameters['rho_d']      #//int(self.Para['Nb'].values)
        self.g = self.parameters['g']      #//int(self.Para['Nb'].values)
        self.L0 = self.parameters['Ls']      #//float(self.Para['Ls'].values)
        self.dt = self.parameters['dtprint']      #//float(self.Para['dt'].values)
        self.g = 1.
        # setting up the data 
        self.tmin = round(tmin,2)
        self.bubbles=list()
        print('reading Pos...')
        self.Pos = pd.read_csv(self.PATH+'Pos.csv');
        # self.Pos = self.Pos[self.Pos.t.ge(self.tmin)]
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
        if self.DIST:
            if 'Dist.csv' in os.listdir(self.PATH):
                file = self.PATH+'Dist.csv'
                print('reading Dist...')
                self.Dist = pd.read_csv(file)
                self.Data = self.Dist[self.Dist.t.ge(self.tmin)]
                self.DataRaw = self.Data.drop(columns=['i','t'])
            if 'Dist_x.csv' in os.listdir(self.PATH) and 'Dist_y.csv' in os.listdir(self.PATH):
                file = self.PATH+'Dist_x.csv'
                print('reading DistX...')
                self.Distx = pd.read_csv(file)
                file = self.PATH+'Dist_y.csv'
                print('reading DistY...')
                self.Disty = pd.read_csv(file)
                self.Dist = np.sqrt(self.Distx**2 + self.Disty**2)
                self.Data = self.Dist[self.Dist.t.ge(self.tmin)]
                self.DataRaw = self.Data.drop(columns=['i','t'])
        
        # samples parameters 
        self.TMAX = self.Vels.t.values[-1]
        self.range = 10.
        self.pair_number = self.Nb*(self.Nb-1)/2
        self.PDF,self.PDFx,self.CDF = {},{},{}
        
        # Norme de U
        self.U = math.sqrt(self.g*self.D*(self.rho_f-self.rho_d)/self.rho_f)
        self.Pos['distmin'] = np.sqrt(self.Pos.distminX.values**2 + self.Pos.distminY.values**2)
        return self 
    
    # @timedeco
    def calcul_all(self):
        print('Calcul_all...')
        if self.DIST:
            if self.RD : self.calcul_radial_distrib(Dr=self.D/10)
            self.calcul_contact_frequency()
            self.calcul_contact_time()
        # self.UprimeUprime_f()
        self.UprimeUprime_p()
        if self.parameters['dimension']==2:
            self.calcul_relative_vel()
            self.calcul_acc()
            # self.calcul_nrj()
            self.calcul_stats_Vels() 
            # self.PFPs()
        elif self.parameters['dimension']==3:
            self.calcul_relative_vel3D()
            self.calcul_acc3D()
            self.calcul_nrj3D()
            self.calcul_stats_Vels3D() 
            
        self.calcul_nearest_stats()
        self.create_bubbles()
        return self
    
    def create_bubbles(self):
        for i in range(self.Nb):
            self.bubbles.append(Bubble(i))
            self.bubbles[i].Pos = self.Pos[self.Pos.tag == i].sort_values(['t']).copy()
        return self
    
    def calcul_nrj(self):
        #energie kinetic
        self.Pos['E_c'] = 1/2*self.Pos['vol']*self.parameters['rho_d']*(self.Pos['vx']**2+self.Pos['vy']**2)
        self.Pos['E_crel'] = 1/2*self.Pos['vol']*self.parameters['rho_d']*(self.Pos['vx_rel']**2+self.Pos['vy_rel']**2)
        #Energy kinetic of rotation
        # self.Pos['E_co'] = 1/2*self.Pos['J_z']*(self.Pos['omegaz']**2)
        # self.Pos['E_s'] = self.parameters['sigma'] * surface area
    
    def calcul_nrj3D(self):
        #energie kinetic
        self.Pos['E_c'] = 1/2*self.Pos['vol']*self.parameters['rho_d']*(self.Pos['vx']**2+self.Pos['vy']**2+self.Pos['vz']**2)
        self.Pos['E_crel'] = 1/2*self.Pos['vol']*self.parameters['rho_d']*(self.Pos['vx_rel']**2+self.Pos['vy_rel']**2+self.Pos['vz_rel']**2)
        self.Pos['E_coz'] = 1/2*self.Pos['J_z']*(self.Pos['omegaz']**2)
        self.Pos['E_coy'] = 1/2*self.Pos['J_y']*(self.Pos['omegay']**2)
        self.Pos['E_cox'] = 1/2*self.Pos['J_x']*(self.Pos['omegax']**2)
        self.Pos['E_co'] = np.sqrt(self.Pos['E_cox']**2+self.Pos['E_coy']**2+self.Pos['E_coz']**2) 
        
        #Energy kinetic of rotation
        # self.Pos['E_s'] = self.parameters['sigma'] * surface area
    

            
    # def calcul_nearest_stats(self,tmin =200,bins=100):
    #     self.N = []
    #     self.DISTMIN = []
    #     self.STDMIN = []
    #     end = self.parameters['DP']
    #     D = self.parameters['D']
    #     data = self.Pos[self.Pos.t.ge(tmin)]
    #     H, valuex = np.histogram(data.values, bins =bins, range=[0,end])
    #     x1 = (valuex[:-1] + valuex[1:])/2
    #     number_of_samples = len(data.index) 
    #     density_of_particles = 1/self.L0**2
    #     weight = math.pi*(valuex[1:]**2-valuex[:-1]**2) *  self.pair_number * number_of_samples * density_of_particles
    #     self.d_nrb = x1
    #     self.Pd_nrb = H
    #     self.Pd_nrbn = H/weight
    #     return self
            
    def calcul_nearest_stats(self,bins= 100,tmin =100):
        self.DISTMIN = []
        self.STDMIN = []
        Dmax = self.parameters['DP']
        Data = self.Pos[self.Pos.t.ge(self.tmin)]
        H, valuex = np.histogram(Data.distmin.values, bins =bins,range=[0,Dmax])
        x1 = (valuex[:-1] + valuex[1:])/2
        number_of_samples = len(Data.index) 
        density_of_particles = 1/self.L0**2
        weight = math.pi*(valuex[1:]**2-valuex[:-1]**2) *  self.pair_number * number_of_samples * density_of_particles
        self.d_nbr = x1
        self.Pd_nbr = H
        self.Pd_nbrn = H/weight
        return self
    
    def calcul_theta_stats(self,tmin =200):
        self.N = []
        self.DISTMIN = []
        self.STDMIN = []
        end = self.parameters['DP']
        D = self.parameters['D']
        data = self.Pos[self.Pos.t.ge(tmin)]
        for C in np.linspace(0.75*D,end,25):
            constrain2  = data.distmin< C + self.parameters['D']/8  
            constrain1  = data.distmin> C - self.parameters['D']/8 
            c  = constrain1 & constrain2
            self.DISTMIN.append(data[c].mean())
            self.STDMIN.append(data[c].std())
            self.N.append(len(data[c]))
        return self
            
      
    def calcul_acc(self):
        self.Pos = self.Pos.sort_values(['tag','t'])
        self.Pos['ax'] = self.Pos['vx'].diff()/self.dt
        self.Pos['ax_rel'] = self.Pos['vx_rel'].diff()/self.dt
        self.Pos['ay'] = self.Pos['vy'].diff()/self.dt
        self.Pos['ay_rel'] = self.Pos['vy_rel'].diff()/self.dt
        self.Pos['a_rel'] = np.sqrt(self.Pos['ax_rel']**2+self.Pos['ay_rel']**2) 
        self.Pos['aomegaz'] = self.Pos['omegaz'].diff()/self.dt
        self.Pos['aomegaz_rel'] = self.Pos['omegaz_rel'].diff()/self.dt
        self.Pos = self.Pos.sort_values(['t','tag'])
        self.Pos = self.Pos[self.Pos.t != 0.1]
        return self
    def calcul_acc_n(self):
        self.Pos = self.Pos.sort_values(['tag','t'])
        self.Pos['an_rel'] =  self.Pos['vn_rel'].diff()/self.dt
        self.Pos['aomegaz'] = self.Pos['omegaz'].diff()/self.dt
        self.Pos['aomegaz_rel'] = self.Pos['omegaz_rel'].diff()/self.dt
        self.Pos = self.Pos.sort_values(['t','tag'])
        self.Pos = self.Pos[self.Pos.t != 0.1]
        
    
    
    def calcul_acc3D(self):
        self.Pos = self.Pos.sort_values(['tag','t'])
        self.Pos['ax'] = self.Pos['vx'].diff()/self.dt
        self.Pos['ax_rel'] = self.Pos['vx_rel'].diff()/self.dt
        self.Pos['ay'] = self.Pos['vy'].diff()/self.dt
        self.Pos['ay_rel'] = self.Pos['vy_rel'].diff()/self.dt
        self.Pos['az'] = self.Pos['vz'].diff()/self.dt
        self.Pos['az_rel'] = self.Pos['vz_rel'].diff()/self.dt
        self.Pos['a_rel'] = np.sqrt(self.Pos['ax_rel']**2+self.Pos['ay_rel']**2+self.Pos['az_rel']**2) 
        self.Pos['aomegaz'] = self.Pos['omegaz'].diff()/self.dt
        self.Pos['aomegay'] = self.Pos['omegay'].diff()/self.dt
        self.Pos['aomegax'] = self.Pos['omegax'].diff()/self.dt
        self.Pos['aomegaz_rel'] = self.Pos['omegaz_rel'].diff()/self.dt
        self.Pos = self.Pos.sort_values(['t','tag'])
        self.Pos = self.Pos[self.Pos.t != 0.1]

    def calcul_relative_vel(self):
        self.Pos = self.Pos.sort_values(['t','tag'])
        dtinv = int(1./self.dt)
        self.Pos['tagminall'] = ( self.Pos.tagmin + (self.Pos.t - self.dt) * dtinv * self.Nb)
        self.Pos = self.Pos.astype({"tagminall": int})
        self.Pos['vx_rel'] = self.Pos.vx.values - self.Pos.iloc[self.Pos.tagminall].vx.values
        self.Pos['vy_rel'] = self.Pos.vy.values - self.Pos.iloc[self.Pos.tagminall].vy.values
        self.Pos['v_rel'] = np.sqrt(self.Pos['vx_rel']**2+self.Pos['vy_rel']**2)
        self.Pos['omegaz_rel'] = self.Pos.omegaz.values - self.Pos.iloc[self.Pos.tagminall].omegaz.values
        self.Pos['nx'] = self.Pos.distminX/self.Pos.distmin
        self.Pos['ny'] = self.Pos.distminY/self.Pos.distmin
        self.Pos['vn_rel'] = (self.Pos.vx.values*self.Pos.nx.values+self.Pos.vy.values*self.Pos.ny.values) - (self.Pos.iloc[self.Pos.tagminall].vx.values*self.Pos.iloc[self.Pos.tagminall].nx.values +self.Pos.iloc[self.Pos.tagminall].vy.values*self.Pos.iloc[self.Pos.tagminall].ny.values)
        # self.Pos['deq'] = np.sqrt(4.*self.Pos.vol.values/math.pi) * np.sqrt(4.*self.Pos.iloc[4.*self.Pos.tagminall].vol.values/math.pi) *2./(np.sqrt(4.*self.Pos.vol.values/math.pi) + np.sqrt(4.*self.Pos.iloc[self.Pos.tagminall].vol.values/math.pi))
        deq = 1.
        self.Pos['We'] = self.parameters['rho_f'] * self.Pos.vn_rel.values**2 /2. /self.parameters['sig']  *deq     
        self.Pos['td'] = self.parameters['rho_f'] * np.abs(self.Pos.v_rel.values) /8. /self.parameters['sig']  *deq**2    / math.sqrt(self.D/self.g)
        self.Pos['beta'] = np.arccos(np.abs((self.Pos.distminX.values*self.Pos.vx_rel.values+self.Pos.distminY.values*self.Pos.vy_rel.values)/(self.Pos.distmin*self.Pos.v_rel)))
        self.Pos['theta'] = np.arccos(np.abs((self.Pos.distminX.values)/(self.Pos.distmin)))
        return self
        

    def calcul_relative_vel3D(self):
        self.Pos = self.Pos.sort_values(['t','tag'])
        dtinv = int(1./self.dt)
        self.Pos['tagminall'] = ( self.Pos.tagmin + (self.Pos.t - self.dt) * dtinv * self.Nb)
        self.Pos = self.Pos.astype({"tagminall": int})
        self.Pos['x_rel'] = self.Pos.x.values - self.Pos.iloc[self.Pos.tagminall].x.values
        self.Pos['y_rel'] = self.Pos.y.values - self.Pos.iloc[self.Pos.tagminall].y.values
        self.Pos['z_rel'] = self.Pos.z.values - self.Pos.iloc[self.Pos.tagminall].z.values
        self.Pos['vx_rel'] = self.Pos.vx.values - self.Pos.iloc[self.Pos.tagminall].vx.values
        self.Pos['vy_rel'] = self.Pos.vy.values - self.Pos.iloc[self.Pos.tagminall].vy.values
        self.Pos['vz_rel'] = self.Pos.vz.values - self.Pos.iloc[self.Pos.tagminall].vz.values
        self.Pos['v_rel'] = np.sqrt(self.Pos['vx_rel']**2+self.Pos['vy_rel']**2+self.Pos['vz_rel']**2)
        self.Pos['omegax_rel'] = self.Pos.omegax.values - self.Pos.iloc[self.Pos.tagminall].omegax.values
        self.Pos['omegay_rel'] = self.Pos.omegay.values - self.Pos.iloc[self.Pos.tagminall].omegay.values
        self.Pos['omegaz_rel'] = self.Pos.omegaz.values - self.Pos.iloc[self.Pos.tagminall].omegaz.values
        self.Pos['omega_rel'] = np.sqrt(self.Pos['omegax_rel']**2+self.Pos['omegay_rel']**2+self.Pos['omegaz_rel']**2)
        self.Pos['distmin'] = np.sqrt(self.Pos.distminX.values**2 + self.Pos.distminY.values**2+ self.Pos.distminZ.values**2)
            
    @timedeco
    def collision_corr(self):
        # data[[index for index in data.keys() if '1' in index.split('b')]]
        return 1

    # def linear_differences(self):
        


    def calcul_radial_distrib(self,bins= 100,tmin =200):
        Dmax = self.D*3
        self.Data = self.Dist[self.Dist.t.ge(self.tmin)]
        self.DataRaw = self.Data.drop(columns=['i','t'])
        H, valuex = np.histogram(self.DataRaw.stack().values, bins =bins,range=[0,Dmax])
        x1 = (valuex[:-1] + valuex[1:])/2
        number_of_samples = len(self.DataRaw.index) 
        density_of_particles = 1/self.L0**2
        weight = math.pi*(valuex[1:]**2-valuex[:-1]**2) *  self.pair_number * number_of_samples * density_of_particles
        self.r = x1
        self.Pr = H/weight
        return self
    
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

    def calcul_contact_frequencyG(self,G=1):
        self.HzG = []
        DataB = self.DataRaw < self.D * G
        Nc=((DataB != DataB.shift(1)) & DataB).sum().sum()# number of contact
        Dt = self.Data.t.values[-1] - self.tmin
        vol = self.L0**2
        self.HzG = Nc /(Dt*self.Nb)  * math.sqrt(self.D/self.parameters['g'])
        return self
    
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
        self.Tcmean = np.mean(self.time_of_contacts)
        # self.PDFx['ct'],self.CDF['ct'],self.PDF['ct'] = self.calcul_PDF(self.time_of_contacts)
        return self

    def calcul_stats_Vels(self,tmin=200):
        self.Vels['dVy'] = self.Vels['vydrops']-self.Vels['vyfluid']
        self.Vels['dVx'] = self.Vels['vxdrops']-self.Vels['vxfluid']
        self.Stats = self.Vels[self.Vels.t.ge(self.tmin)]
        self.dVy = self.Vels.dVy[self.Vels.t.ge(tmin)].describe()
        self.F =self.parameters['g'] * (self.parameters['rho_f'] - self.parameters['rho_d']) *self.parameters['PHI']*self.parameters['Ls']**2/self.Nb
        self.F_mu =  self.F / (self.parameters['mu_f']*self.dVy['mean']*self.parameters['D'])
        self.F_mustd =  self.dVy['std']*self.F / (self.parameters['mu_f']*self.dVy['mean']**2*self.parameters['D'])
        self.F_rho =  self.F / (0.5*self.parameters['rho_f']*self.dVy['mean']**2*self.parameters['D'])
        self.F_rhostd =  self.dVy['std']*self.F / (0.5*self.parameters['rho_f']*self.dVy['mean']**3*self.parameters['D'])
        return self
        
    def calcul_stats_Vels3D(self):
        self.Vels['dVy'] = self.Vels['vydrops']-self.Vels['vyfluid']
        self.Vels['dVx'] = self.Vels['vxdrops']-self.Vels['vxfluid']
        self.Vels['dVz'] = self.Vels['vzdrops']-self.Vels['vzfluid']
        self.Stats = self.Vels[self.Vels.t.ge(self.tmin)]
        
        # valocities fluctuation 
     
    
    def calcul_PDF(self,list,range = None):
        bins = 100
        H,X1=np.histogram(list,bins=bins,range=range)
        dx = X1[1]-X1[0]
        CDF = np.cumsum(H)/len(list)
        PDFx,PDF = X1[1:],np.gradient(CDF,dx)
        return PDFx,CDF,PDF

    def UprimeUprime_p(self,tmin=100):
        """calculation of the particular average pseudo turbulent tensor
        """
        # can be done with describe module of dataframe objject 
        self.Pos['Uprime_x'] = self.Pos['vx'] - self.Pos['vx'].mean()
        self.Pos['Uprime_y'] = self.Pos['vy'] - self.Pos['vy'].mean()
        self.Pos['UU_pxx'] = self.Pos['Uprime_x']**2
        self.Pos['UU_pyy'] = self.Pos['Uprime_y']**2
        self.Pos['UU_pxy'] = self.Pos['Uprime_x']*self.Pos['Uprime_y']
        self.Pos['Omegap'] = self.Pos['omegaz'] - self.Pos['omegaz'].mean()
        self.Pos['OU_pzx'] = self.Pos['Omegap']*self.Pos['Uprime_x']
        self.Pos['OU_pzy'] = self.Pos['Omegap']*self.Pos['Uprime_y']
        self.UU_pxx = self.Pos[self.Pos.t.ge(tmin)].UU_pxx.describe()
        self.UU_pyy = self.Pos[self.Pos.t.ge(tmin)].UU_pyy.describe()
        self.UU_pxy = self.Pos[self.Pos.t.ge(tmin)].UU_pxy.describe()
        self.UU_fxx = self.Vels[self.Vels.t.ge(tmin)].xp.describe()
        self.UU_fyy = self.Vels[self.Vels.t.ge(tmin)].yp.describe()
        self.OU_pzx = self.Pos[self.Pos.t.ge(tmin)].OU_pzx.describe()
        self.OU_pzy = self.Pos[self.Pos.t.ge(tmin)].OU_pzy.describe()
        # normaliser par U**2 
        return self
        

    def PFPs(self):
        x = ['x','y']
        X = ['X','Y']
        for i in X:
            for j in x:
                self.Pos['r'+i+'f'+j] = self.Pos['distmin'+i] * self.Pos['a'+j] * self.Pos['vol'] * self.parameters['rho_d']
    
     