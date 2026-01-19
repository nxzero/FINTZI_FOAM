import pandas as pd
import better_py_functions.PPcyl as pc
import better_py_functions.PPsus as ps
from py_functionsBA.header import loadbar
import os
import numpy as np
RESULT_DIR='/home/irsrvshare1/R04/RheoPipe/fintzin/CYLINDERS_PROJECT.backup/results/'
class PPcla():
    def __init__(self,namestud='',forcedir='forces',resultdir = RESULT_DIR):
    # def Postprocess(self, namedir):
        self.namestud = namestud
        self.resultsdir = resultdir
        if not namestud == '':
            self.namedir = self.resultsdir+namestud+'/'
            self.forcedir = forcedir
            self.names = pd.Series(os.listdir(self.namedir))
            self.names = self.names[self.names.str.replace('.','').str.isdigit()]#
            self.names = self.names.astype(float).sort_values(ascending=True).astype(str)
            for i,n in enumerate(self.names) :
                if n[-2:] =='.0':
                    n = n[:-2]
                self.names.iloc[i] = n
            self.parameters_list = [float(name) for name in self.names]
        else:
            self.names = []
            self.initlist()
        self.momentsAD = pd.DataFrame()
        self.momentsD = pd.DataFrame()
        self.moments = pd.DataFrame()
        self.forcesAD = pd.DataFrame()
        self.forcesD = pd.DataFrame()
        self.forces = pd.DataFrame()
        self.firstmomentsD = pd.DataFrame()
        self.firstmoments = pd.DataFrame()
        self.Turb = pd.DataFrame()
        self.parameters = pd.DataFrame()
        
    def load_PPcyl(self):
        self.Cases = list()
        for i,name in enumerate(self.names):
            a = pc.PPcyl(self.namestud+'/'+name,forces_dir=self.forcedir,resultdir=self.resultsdir)
            a.get_forces()
            a.get_moments()
            self.Cases.append(a)
            
        self.make_list_para()
        self.make_list_forces()
        self.make_list_moments()
        print(self.namestud,":",self.parameters_list)
        return self

    def make_list_forces(self):
        for c in  self.Cases:
            # self.forcesAD = self.forcesAD.append(c.forcesAD.iloc[-1:], ignore_index = True)
            # self.forcesD = self.forcesD.append(c.forcesD.iloc[-1:], ignore_index = True)
            # self.forces = self.forces.append(c.forces.iloc[-1:], ignore_index = True)
            self.forcesAD = pd.concat([self.forcesAD, c.forcesAD.iloc[-1:]], ignore_index=True)
            self.forcesD = pd.concat([self.forcesD, c.forcesD.iloc[-1:]], ignore_index=True)
            self.forces = pd.concat([self.forces, c.forces.iloc[-1:]], ignore_index=True)

    def make_list_moments(self):
        for c in  self.Cases:
            # self.momentsAD = self.momentsAD.append(c.momentsAD.iloc[-1:], ignore_index = True)
            # self.momentsD = self.momentsD.append(c.momentsD.iloc[-1:], ignore_index = True)
            # self.moments = self.moments.append(c.moments.iloc[-1:], ignore_index = True)
            self.momentsAD = pd.concat([self.momentsAD, c.momentsAD.iloc[-1:]], ignore_index=True)
            self.momentsD = pd.concat([self.momentsD, c.momentsD.iloc[-1:]], ignore_index=True)
            self.moments = pd.concat([self.moments, c.moments.iloc[-1:]], ignore_index=True)


    def make_list_firstmoments(self):
        for c in  self.Cases:
            # self.firstmomentsD = self.firstmomentsD.append(c.firstmomentsD.iloc[-1:], ignore_index = True)
            # self.firstmoments = self.firstmoments.append(c.firstmoments.iloc[-1:], ignore_index = True)
            self.firstmoments = pd.concat([self.firstmoments, c.firstmoments.iloc[-1:]], ignore_index=True)
            self.firstmomentsD = pd.concat([self.firstmomentsD, c.firstmomentsD.iloc[-1:]], ignore_index=True)


    def load_PPsus(self):
        self.Cases = list()
        for i,name in enumerate(self.names):
            print(self.namestud+name)
            a = ps.PPsus(self.namestud+'/'+name,forces_dir=self.forcedir,resultdir=self.resultsdir)
            a.load_PPcyl()
            a.calcDp()
            a.get_Turb()
            self.Cases.append(a)
            loadbar(i,len(self.names),prefix='READING_DATAS')
        self.make_list_para()
        self.make_list_Dp()
        self.make_list_UU()
        self.make_avg_list_forces()
        self.make_avg_list_firstmoments()
        return self
        
        
    def make_list_para(self):
        for c in  self.Cases:
            para = pd.DataFrame([c.parameters])
            # self.parameters = self.parameters.append(para)
            self.parameters = pd.concat([self.parameters, para], ignore_index=True)
            

    def make_list_UU(self):
        for c in  self.Cases:
            # self.Turb = self.Turb.append(c.Turb.iloc[-1:], ignore_index = True)
            self.Turb = pd.concat([self.Turb, c.Turb.iloc[-1:]], ignore_index=True)
            
    def make_list_Dp(self):
        self.Dp = list()
        self.Re_eq = list()
        for c in  self.Cases:
            self.Dp.append(c.Dp)
            self.Re_eq.append(c.Re_eq)
        self.Dp = pd.Series(self.Dp)
        self.Re_eq = np.array(self.Re_eq)

    def make_avg_list_forces(self):
        for c in  self.Cases:
            # self.forcesD = self.forcesD.append(c.forcesD.mean(), ignore_index = True)
            # self.forces = self.forces.append(c.forces.mean(), ignore_index = True)
            self.forcesD = pd.concat([self.forcesD, c.forcesD.mean()], ignore_index=True)
            self.forces = pd.concat([self.forces, c.forces.mean()], ignore_index=True)
            

    def make_avg_list_firstmoments(self):
        for c in  self.Cases:
            # self.firstmomentsD = self.firstmomentsD.append(c.firstmomentsD.mean(), ignore_index = True)
            # self.firstmoments = self.firstmoments.append(c.firstmoments.mean(), ignore_index = True)
            self.firstmomentsD = pd.concat([self.firstmomentsD, c.firstmomentsD.mean()], ignore_index=True)
            self.firstmoments = pd.concat([self.firstmoments, c.firstmoments.mean()], ignore_index=True)
            
# from better_py_functions.PPcla import *
# a = Study('TP_clean/200C/CHI_5/PHI_25bis/',resultdir='../../CYLINDERS_PROJECT/results/')
# a.load_PPsus()