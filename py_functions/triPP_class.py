from py_functions.PP_class import *
import numpy as np
from numpy import linalg as LA
import math
import sys
class triStudy(Study):
    def __init__(self, namestud):
        Study.__init__(self, namestud)
        self.noc = self.parameters['noc']
        self.namestud = namestud
        self.Cyls = {}
        self.SUM ={}
        self.DIST ={}
        self.DISTC ={}
        self.AVG ={}
        self.DEV ={}
        self.frequency = 10
        self.DISTC = {}
        self.Turb = {}
        self.Dp = {}
        self.OT = {}
        self.cylsDeleted = []
    def get_studies_last_step(self,name):
        self.dircase = self.namedir+str(name)
        
        sys.path.append('./'+self.dircase+'/')
        import Cyls_for_this_study as cylsStud
        from importlib import reload
        reload(cylsStud)
        self.cylsStud = cylsStud
        sys.path.remove('./'+self.dircase+'/')

        self.Cyls[name] = []
        self.noc = int(self.noc)
        if self.parameters['para_name'] == 'noc':
            self.noc = int(name)

        pp = np.zeros((3,3))
        for i in range(self.noc):
            maindir = 'forces_' + str(i)
            self.dirs = [di for di in os.listdir(self.dircase) if str(i) in di.split('_') and 'forces' in di.split('_')]
            if self.dirs == []:
                self.cylsDeleted.append(i)
                continue
            Id = '_'.join(self.dirs[0].split('_')[1:])

            self.Cyls[name].append(Study(self.namestud))
            self.Cyls[name][-1].name = name
            self.Cyls[name][-1].get_study_last_step(name,self.dirs[0])
            for other in self.dirs[1:]:
                cyl = Study(self.namestud)
                cyl.get_study_last_step(name,other)
                self.Cyls[name][-1] += cyl
            Index = cylsStud.Ids.index(Id)
            ori = np.array(cylsStud.ori[Index])
            ID = np.array(cylsStud.Ids[Index])
            pp += np.outer(ori,ori)
            self.Cyls[name][-1].ori = ori
            self.Cyls[name][-1].ID = ID
            # fix the real angle of the cylinder
            self.Cyls[name][-1].studies_datas[name]['parameters']['Theta'] = np.arccos(abs(ori[0])) /(math.pi *2)*360 
            # fix the axis on which the torque will be calculated 
            self.Cyls[name][-1].calcul(name)
        self.OT[name] = pp/self.noc



    def get_forces_and_troques(self):
        print(self.namestud,self.para_name,self.names)
        for name in self.names:
            self.get_studies_last_step(name)
            self.calcTurb(name)
            self.calcDp(name)
        for name in self.names:
            for cyl in self.Cyls[name]:
                cyl.initlist()
                cyl.makelist(name)
            self.calcul_sum(name)
            self.distribution_theta(name)
            self.average_theta(name)
            self.standard_deviation_theta(name)

    def calcul_sum(self,name):
        self.SUM[name] = Study(self.namestud)
        self.SUM[name].name = name
        self.SUM[name].get_study_last_step(name,self.dirs[0])
        self.SUM[name].set_to_0(self.Cyls[name][0])
        for cyl in self.Cyls[name]:
            self.SUM[name].addall(cyl)
        self.SUM[name].initlist()
        self.SUM[name].makelist(name)

        self.SUM[name].AVG = Study(self.namestud)
        self.SUM[name].AVG.name = name
        self.SUM[name].AVG.get_study_last_step(name,self.dirs[0])
        self.SUM[name].AVG.set_to_0(self.Cyls[name][0])
        for cyl in self.Cyls[name]:
            self.SUM[name].AVG.addall(cyl)
        self.SUM[name].AVG.divall(float(self.noc))
        
        self.SUM[name].AVG.initlist()
        self.SUM[name].AVG.makelist(name)

    def distribution_theta(self,name):
        self.DIST[name] = {}
        self.DIST[name]['ranges'] = np.linspace(0,90,self.frequency)
        self.DIST[name]['Theta'] = np.zeros(len(self.DIST[name]['ranges'])-1)
        i = 0
        for theta_min,theta_max in zip(self.DIST[name]['ranges'][:-1],self.DIST[name]['ranges'][1:]):
            for cyl in self.Cyls[name]:
                theta = cyl.studies_datas[name]['parameters']['Theta']
                if theta <= theta_max and theta > theta_min:
                    self.DIST[name]['Theta'][i] += 1
            i += 1
        # self.DISTC[name] = self.DIST[name].copy()
        # i = 0
        # for theta_min,theta_max in zip(self.DIST[name]['ranges'][:-1],self.DIST[name]['ranges'][1:]):
        #     for j in range(i):
        #         self.DISTC[name]['Theta'][i] += self.DIST[name]['Theta'][j]
        #     i += 1


    def average_theta(self,name):
        self.AVG[name] = {}
        self.AVG[name]['ranges'] = np.linspace(0,90,self.frequency)
        self.AVG[name]['Theta'] = list()
        i = 0
        for theta_min,theta_max in zip(self.AVG[name]['ranges'][:-1],self.AVG[name]['ranges'][1:]):
            self.AVG[name]['Theta'].append(Study(self.namestud))
            self.AVG[name]['Theta'][i].name = name
            self.AVG[name]['Theta'][i].get_study_last_step(name,self.dirs[0])
            self.AVG[name]['Theta'][i].set_to_0(self.Cyls[name][0])
            for cyl in self.Cyls[name]:
                theta = cyl.studies_datas[name]['parameters']['Theta']
                if theta <= theta_max and theta > theta_min:
                    self.AVG[name]['Theta'][i].addall(cyl)
            self.AVG[name]['Theta'][i].divall(self.DIST[name]['Theta'][i])
            self.AVG[name]['Theta'][i].initlist()
            self.AVG[name]['Theta'][i].makelist(name)
            i += 1
        
        self.AVG[name]['ranges'] = self.AVG[name]['ranges'][:-1] + abs(theta_min-theta_max)/2


        #################### add total avg 
    def standard_deviation_theta(self,name):
        self.DEV[name] = {}
        self.DEV[name]['ranges'] = np.linspace(0,90,self.frequency)
        self.DEV[name]['Theta'] = list()
        i = 0
        for theta_min,theta_max in zip(self.DEV[name]['ranges'][:-1],self.DEV[name]['ranges'][1:]):
            self.DEV[name]['Theta'].append(Study(self.namestud))
            self.DEV[name]['Theta'][i].name = name
            self.DEV[name]['Theta'][i].get_study_last_step(name,self.dirs[0])
            self.DEV[name]['Theta'][i].set_to_0(self.Cyls[name][0])
            for cyl in self.Cyls[name]:
                theta = cyl.studies_datas[name]['parameters']['Theta']
                if theta <= theta_max and theta > theta_min:
                    self.DEV[name]['Theta'][i].addall((cyl.suball(self.AVG[name]['Theta'][i]))**2)
            self.DEV[name]['Theta'][i].sqrt()
            self.DEV[name]['Theta'][i].divall(np.sqrt(self.DIST[name]['Theta'][i]))
            self.DEV[name]['Theta'][i].studies_datas[name]['parameters']['Theta'] = (theta_min+theta_max)/2
            self.DEV[name]['Theta'][i].initlist()
            self.DEV[name]['Theta'][i].makelist(name)
            i += 1
        
        self.DEV[name]['ranges'] = self.DEV[name]['ranges'][:-1] + abs(theta_min-theta_max)/2
        #################### add total DEV
                
    def calcTurb(self,name):
        #definition des parametres
        dirF = self.namedir+str(name)+'/Turb1/'
        dirs = os.listdir(dirF)
        dirs = filter(lambda x : x!='' and x !='.' and x!='0.' , [re.sub("[^0123456789\.]","",dir) for dir in dirs])
        dirs =[dire for dire in dirs if dire in os.listdir(dirF)]
        dirs.sort(key = lambda x : float(x))#range dans l'ordre croissant 
        lastdir = dirs[-1]
        
        Forces_path = dirF + lastdir +'/Turb1.dat'
        Forces = open(Forces_path,'r').readlines()[-1].split()
        Forces = [float(re.sub("[^0123456789\.e+-]","",force)) for force in Forces]
        self.Turb[name] = {
            'XX' : Forces[1],
            'YY' : Forces[2],
            'ZZ' : Forces[3],
            'XY' : Forces[4],
            'XZ' : Forces[5],
            'ZY' : Forces[6]
        }
    def calcDp(self,name):
        #definition des parametres
        dirF = self.namedir+str(name)+'/system/Dp'
        Dp = open(dirF,'r').readlines()[0].split()[1]
        Dp = float(re.sub("[^0123456789\.e+-]","",Dp)) 
        self.Dp[name] = Dp * self.parameters['Rho']

