import re
import os
import math
import sys
import os
import re
import numpy as np
import numpy.linalg as LA
import pandas as pd
RESULT_DIR='/home/irsrvshare1/R04/RheoPipe/fintzin/CYLINDERS_PROJECT.backup/results/'
class PPcyl():
    def __init__(self,namecase,forces_dir='forces',resultdir = RESULT_DIR,para = dict()):
    # def Postprocess(self, namedir):
        self.namecase = namecase
        self.resultsdir = resultdir
        self.namedir = self.resultsdir+namecase+'/'
        # print(self.namedir)
        sys.path.append(self.namedir)
        import parameters_for_this_study as pstud
        from importlib import reload
        reload(pstud)
        pstud.parameters.update(para)
        self.parameters = pstud.parameters
        sys.path.remove(self.namedir)
        self.forces_dir = self.resultsdir + self.namecase +'/'+forces_dir +'/'
        self.para_name = self.parameters['para_name']
        self.ori = np.array([6,6,6])
        ##### definition of the scalling para

        # to get drag coefficients 
        Ukey = 'U'
        if self.parameters['U'] == 0:
            Ukey = 'Uz'
        if 'omega' not in self.parameters.keys():
            self.parameters['omega'] = 0
        self.geometrycalparameters = 1./2.*self.parameters['Rho']*self.parameters['r']*2*self.parameters[Ukey]**2*self.parameters['length']
        self.F_h = 6*self.parameters['eta']*math.pi*self.parameters[Ukey]*(3./4.*self.parameters['r']**2*self.parameters['length'])**(1./3.)
        # self.F_h = self.parameters['eta']*self.parameters[Ukey]*((self.parameters['r']*2)**3)
        if self.parameters['whichMesh'] in [4,3]:
            omega   = 'omega'
        else:
            omega = 'omega_para'
        try:self.parameters[omega]
        except: omega = 'U'
            
        # self.T_h = 8*self.parameters['eta']*math.pi*self.parameters[omega]*(3./4.*self.parameters['r']**2*self.parameters['length'])
        self.T_h = self.parameters['eta']*self.parameters[omega]*(2*self.parameters['r'])**2*self.parameters['length']
        self.adimensionneur = self.parameters['eta'] * self.parameters[Ukey] *self.parameters['length'] 
        self.adimensionneur2 = self.parameters['eta'] * self.parameters[Ukey] *self.parameters['length']**2 
        ####################################################
        self.Id = ''
        self.pos = np.array([6,6,6]) #random values
        self.ori = np.array([6,6,6]) #random values
        self.forces = pd.DataFrame()
        self.moments = pd.DataFrame()
        self.firstmoments = pd.DataFrame()
        
    def get_forces(self):
        steps = pd.Series(os.listdir(self.forces_dir))
        laststep = max([float(a) for a in steps[steps.str.replace('.','').str.isdigit()]])
        if str(laststep)[-2:] == '.0':
            laststep = str(laststep)[:-2]
        file =self.forces_dir+str(laststep)+'/force.dat'
        self.forces = pd.read_table(file,sep = '\s+',skiprows=[0,1,2,3],header=None)
        self.forces.columns = ['t','x','y','z','Px','Py','Pz','Vx','Vy','Vz']
        for col in self.forces.columns:
            if self.forces[col].dtype == object: 
                self.forces[col] = self.forces[col].str.replace('[()]','',regex=True).astype(float) 
        # drag on a sphere of same vol
        self.forcesAD = self.forces/self.F_h
        self.forcesD = self.forces/self.adimensionneur
    
    def get_moments(self):
        steps = pd.Series(os.listdir(self.forces_dir))
        laststep = max([float(a) for a in steps[steps.str.replace('.','').str.isdigit()]])
        if str(laststep)[-2:] == '.0':
            laststep = str(laststep)[:-2]
        self.moments = pd.read_csv(self.forces_dir+str(laststep)+'/moment.dat',sep = '\s+',skiprows=[0,1,2,3],header=None)
        self.moments.columns = ['t','x','y','z','Px','Py','Pz','Vx','Vy','Vz']
        for col in self.moments.columns:
            if self.moments[col].dtype == object: 
                self.moments[col] = self.moments[col].str.replace('[()]','',regex=True).astype(float)

        self.momentsD = self.moments/self.adimensionneur2
        self.momentsAD = self.moments/self.T_h
        
    def get_firstmoments(self):
        steps = pd.Series(os.listdir(self.forces_dir))
        laststep = max([float(a) for a in steps[steps.str.replace('.','').str.isdigit()]])
        if str(laststep)[-2:] == '.0':
            laststep = str(laststep)[:-2]
        self.firstmoments = pd.read_csv(self.forces_dir+str(laststep)+'/firstmoment.dat',sep = '\s+',skiprows=[0,1,2,3],header=None)
        self.firstmoments.columns = ['t','xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz','Pxx','Pxy','Pxz','Pyx','Pyy','Pyz','Pzx','Pzy','Pzz','Vxx','Vxy','Vxz','Vyx','Vyy','Vyz','Vzx','Vzy','Vzz']
        for col in self.firstmoments.columns:
            if self.firstmoments[col].dtype == object: 
                self.firstmoments[col] = self.firstmoments[col].str.replace('[()]','',regex=True).astype(float)
                
        self.firstmomentsD = self.firstmoments/self.adimensionneur2


        
    def __add__(self,other):
        self.forces += other.forces
        self.forcesD += other.forcesD
        self.forcesAD += other.forcesAD
        if not self.moments.empty:
            self.moments += other.moments
            self.momentsAD += other.momentsAD
            self.momentsD += other.momentsD
        if not self.firstmoments.empty:
            self.firstmoments += other.firstmoments
            self.firstmomentsD += other.firstmomentsD
        return self
        
    def __sub__(self,other):
        self.forces -= other.forces
        self.forcesD -= other.forcesD
        self.forcesAD -= other.forcesAD
        if not self.moments.empty:
            self.moments -= other.moments
            self.momentsAD -= other.momentsAD
            self.momentsD -= other.momentsD
        if not self.firstmoments.empty:
            self.firstmoments -= other.firstmoments
            self.firstmomentsD -= other.firstmomentsD
        return self
        
    def tri_perio_calc(self):
        # steps = pd.Series(os.listdir(self.forces_dir))
        # laststep = max([float(a) for a in steps[steps.str.replace('.','').str.isdigit()]])
        # if str(laststep)[-2:] == '.0':
        #     laststep = str(laststep)[:-2]
        # self.firstmoments = pd.read_csv(self.forces_dir+str(laststep)+'/firstmoment.dat',sep = '\s+',skiprows=[0,1,2,3],header=None)
        # self.firstmoments.columns = ['t','xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz','Pxx','Pxy','Pxz','Pyx','Pyy','Pyz','Pzx','Pzy','Pzz','Vxx','Vxy','Vxz','Vyx','Vyy','Vyz','Vzx','Vzy','Vzz']
        # for col in self.firstmoments.columns:
        #     if self.firstmoments[col].dtype == object: 
        #         self.firstmoments[col] = self.firstmoments[col].str.replace('[()]','').astype(float)
        
        
        self.axis_of_rot = np.cross(np.array([1,0,0]), self.ori)
        self.axis_of_rot /= LA.norm(self.axis_of_rot,2)
        self.parameters['thetaRed'] = np.arccos(abs(self.ori[0]))
        self.parameters['theta'] = np.arccos(abs(self.ori[0])) / (math.pi)*180

    def read_feild(self):
        D = self.parameters['r']*2
        L = self.parameters['length']
        chi = self.parameters['ksi']
        xCenter =self.parameters['xCenter']
        yCenter =self.parameters['yCenter']
        zCenter =self.parameters['zCenter']
        eta =self.parameters['eta']
        U =self.parameters['U']
        dH  = D/4. 
        ThetaU = self.parameters['ThetaU']
        Theta = self.parameters['Theta']
        ct = np.cos(ThetaU*math.pi/180)
        ctt = np.cos(Theta*math.pi/180)
        st = np.sin(ThetaU*math.pi/180)
        stt = np.sin(Theta*math.pi/180)
        self.data = pd.read_csv(self.resultsdir+self.namecase+'datas.csv')
        # self.data2 = pd.read_csv(self.resultsdir+self.namecase+'datas2.csv')
        self.s = []
        self.fd = []
        self.fl = []
        self.mz = []
        
        dataB =  self.data[self.data['Block Name'] == 'air_to_spherePlaneFacesbot']
        dataT =  self.data[self.data['Block Name'] == 'air_to_spherePlaneFacestop']
        dataS =  self.data[self.data['Block Name'] == 'air_to_sphere']
        
        
        
        d = dataB
        FD = d['forces:force_0'].sum() * ct + d['forces:force_1'].sum() * st
        FL = d['forces:force_0'].sum() * (-st) + d['forces:force_1'].sum() * ct
        FM = d['forces:moment_2'].sum()
        self.s.append(-L/2)
        self.fd.append(FD)
        self.fl.append(FL)
        self.mz.append(FM)
        N = int(4*chi)
        # listt = list(np.linspace(-L/2+dH/2,L/2-dH/2, 100)
        # listt = [i for j, i in enumerate(listt) if j % 2 == 0]
        # for c in listt:
        dH = D/4.1
        for c in list(np.linspace(-L/2+dH/2,L/2-dH/2 , int(4*chi))):
            print(-L/2+dH/2,L/2-dH/2,c)
            cond1 = ((dataS['C_0']-xCenter)*ctt+(dataS['C_1'] - yCenter)*stt ) <= (c + dH/2.)
            cond2 = ((dataS['C_0']-xCenter)*ctt+(dataS['C_1'] - yCenter)*stt ) >= (c - dH/2.)
            cond = cond1 & cond2
            d = dataS[cond]
            FD = d['forces:force_0'].sum() * ct     + d['forces:force_1'].sum() * st
            FL = d['forces:force_0'].sum() * (-st)  + d['forces:force_1'].sum() * ct
            FM = d['forces:moment_2'].sum()
            self.s.append(c)
            self.fd.append(FD)
            self.fl.append(FL)
            self.mz.append(FM)
        
            
        d = dataT
        FD = d['forces:force_0'].sum() * ct + d['forces:force_1'].sum() * st
        FL = d['forces:force_0'].sum() * (-st) + d['forces:force_1'].sum() * ct
        FM = d['forces:moment_2'].sum()
        self.s.append(L/2.)
        self.fd.append(FD)
        self.fl.append(FL)
        self.mz.append(FM)
        
        self.s = [i /(L) for i in self.s ]
        self.fd = [i /(eta*U*dH) for i in self.fd ]
        self.fl = [i /(eta*U*dH) for i in self.fl ]
        self.mz = [i /(eta*U*dH*L) for i in self.mz ]
        
