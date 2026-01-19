import re
import os
import math
import sys
import os
import re
import numpy as np
import pandas as pd
import better_py_functions.PPcyl as pc
RESULT_DIR='/home/irsrvshare1/R04/RheoPipe/fintzin/CYLINDERS_PROJECT.backup/results/'

class PPsus():
    def __init__(self,namecase,forces_dir='forces',resultdir = RESULT_DIR,para = dict()):
    # def Postprocess(self, namedir):
        self.namecase = namecase
        self.resultsdir = resultdir
        self.namedir = self.resultsdir+namecase+'/'
        sys.path.append('./'+self.namedir)
        import parameters_for_this_study as pstud
        import Cyls_for_this_study as cylsp
        from importlib import reload
        reload(pstud)
        reload(cylsp)
        pstud.parameters.update(para)
        self.Ids = cylsp.Ids
        self.pos = cylsp.pos
        self.ori = cylsp.ori
        self.dim = cylsp.dim # A changer :!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.parameters = pstud.parameters
        sys.path.remove('./'+self.namedir)
        self.forces_dir = forces_dir 
        self.para_name = self.parameters['para_name']
        
        self.forces = pd.DataFrame()  
        self.forcesD = pd.DataFrame()  
        self.moments = pd.DataFrame()  
        self.momentsD = pd.DataFrame()  
        self.firstmoments = pd.DataFrame()  
        self.firstmomentsD = pd.DataFrame()  
        self.adimensionneur = self.parameters['eta'] * self.parameters['U'] *self.parameters['length'] / 2
        self.adimensionneur2 = self.parameters['eta'] * self.parameters['U'] *self.parameters['length']**2 / 4.

    def load_PPcyl(self):
        self.AllCyls = list()
        self.Cyls = list()
        realId =[Id.replace(self.forces_dir+'_','') for Id in os.listdir(self.namedir)]
        for ID,pos,ori,dim in zip(self.Ids,self.pos,self.ori,self.dim):
            if ID in realId:
                force_dir = self.forces_dir+'_'+ID
                #re use PPcyl class for each cyl of the thing 
                #but it takes a lot of time
                a = pc.PPcyl(self.namecase,forces_dir=force_dir,resultdir=self.resultsdir,para={'r':dim[0],'length':dim[1]})
                a.pos  = pos
                a.Id   = ID
                a.ori  = ori
                    
                a.get_forces()
                a.get_firstmoments()
                a.tri_perio_calc()
                self.AllCyls.append(a)
        #adding the differents periodic cylinders together (since thier represent one cyl)
        for i in range(self.parameters['noc']+1):
            list_of_cyls = [c for c in self.AllCyls if c.Id.split('_')[0] == str(i)]
            self.Cyls.append(list_of_cyls[0])
            for c in list_of_cyls[1:]:
                self.Cyls[i] = self.Cyls[i] + c
            self.Cyls[i].forces.t /= len(list_of_cyls)
            self.Cyls[i].firstmoments.t /= len(list_of_cyls)
        self.make_list_cyls()
                
    def make_list_cyls(self):
        para = list()
        para.append(self.Cyls[0].parameters['theta'])
        for c in  self.Cyls:
            if type(self.parameters['wrinteEachTimeStep']) is int: 
                #appending all last values 
                self.forces = self.forces.append(c.forces.iloc[-1:], ignore_index = True)
                self.firstmoments = self.firstmoments.append(c.firstmoments.iloc[-1:], ignore_index = True)
            else:
                #appending mean values in case the case is unsteady
                #i.e. if the simulation ran with pimpleFoam
                self.forces = self.forces.append(c.forces.mean(), ignore_index = True)
                self.firstmoments = self.firstmoments.append(c.firstmoments.mean(), ignore_index = True)
            para.append(c.parameters['theta'])
        self.para = pd.DataFrame({'theta':para}) 
        self.firstmomentsD = self.firstmoments/self.adimensionneur2
        self.forcesD = self.forces/self.adimensionneur
        
     
    def make_list_time(self):
        self.t =  self.Cyls[0].forces.t
        self.forces = self.Cyls[0].forces
        self.firstmoments = self.Cyls[0].firstmoments
        for c in  self.Cyls[1:]:
            self.forces = self.forces + c.forces
            self.firstmoments = self.firstmoments + c.firstmoments
        self.forces = self.forces/len(self.Cyls)
        self.firstmoments = self.firstmoments/len(self.Cyls)
        self.firstmomentsD = self.firstmoments/self.adimensionneur2
        self.forcesD = self.forces/self.adimensionneur
        
    def calcDp(self):
        dirF = self.namedir+'/system/Dp'
        Dp = open(dirF,'r').readlines()[0].split()[1]
        Dp = float(re.sub("[^0123456789\.e+-]","",Dp)) 
        #the length of the domain is x2 so Dp is the pressure gradient 
        self.Dp = Dp * self.parameters['Rho'] / self.parameters['x2']
        
        self.Re_eq = self.parameters['Re'] 

    def get_Turb(self):
        steps = pd.Series(os.listdir(self.namedir+'Turb1/'))
        laststep = max([float(a) for a in steps[steps.str.replace('.','').str.isdigit()]])
        if str(laststep)[-2:] == '.0':
            laststep = str(laststep)[:-2]
        file =self.namedir+'Turb1/'+str(laststep)+'/Turb1.dat'
        self.Turb = pd.read_table(file,sep = '\s+',skiprows=[0,1],header=None)
        self.Turb.columns = ['t','x','y','z','xy','xz','zy']
        # adim fluctuation tensor
        self.Turb = self.Turb/self.parameters['U']**2 
    def get_PP(self):
        T = np.zeros([3,3])
        for c in self.Cyls:
            T += np.outer(c.ori,c.ori)
        self.T = T / len(self.Cyls)
# PHI20_chi5 = PPsus('TP_clean/200C/CHI_5/PHI_20bis/80',resultdir='../../CYLINDERS_PROJECT/results/')
# PHI20_chi5.load_PPcyl()
# PHI20_chi5.get_PP()
# print(PHI20_chi5.T)

