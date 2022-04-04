import numpy as np
import sys
import math 
import os
import filecmp
import shutil
class AP():
    def __init__(self,para_name: str,parameters_list: list,namestud: str):
        self.para_name = para_name
        self.parameters_list = parameters_list
        self.resultsdir = "../results/"
        self.namestud=namestud
        self.namedir 	= self.resultsdir+namestud.strip()+'/'
        self.name_of_C_file = 'rising-suspension'
        self.MPI = True
        self.nProc = 9
        ## Defaults adimensional parameters 
        # physical parameters
        self.Bo  = 0.5              # Bond number
        self.Ga = 55                # Galileo number
        self.rho_r = 1.1668599      # ratio of rho
        self.mu_r = 0.4237805       # ratio between mu
        self.PHI = 0.15             # concentration phi
        self.N = 10                 # sqrt of the number of bubbles 
        # numerical parameters 
        self.ndcmin         = 40
        self.TMAX           = 2000         # end time 
        self.dtprint        = 0.1          # writing time step
        self.FIRSTSTEP      = self.dtprint * 1
        self.SAVESTEP       = self.TMAX/3
        self.MOVIES         = 1
        # geometrical parameters 
        self.dimension = 2                # set the dimension of the simulation 
        self.nopara    = ''
        
        ## constant dimensional parameters 
        self.g = 1.                 # gravity 
        self.r1 = 1.
        self.D = 1.                 # diameter of bubbles
        self.rho_d = 1.             # density of the dispersed phase 
        
        ## depending parameters 
        # Physiical parameters
        if self.dimension == 2:                                                   # if it is in 2d
            self.Nb = self.N**2                              # number of cells 
            self.Ls = self.N*self.D*math.sqrt(math.pi/(4.*self.PHI))        # size of the domain
        if self.dimension == 3:                                                   # if it is in 3d
            self.Nb = self.N**3                         # number of cells 
            self.Ls = self.N*self.D*math.pow(math.pi/(6.*self.PHI),1./3.)   # size of the domain
        self.mu_d = math.sqrt(abs(1-self.rho_r)*self.rho_r*self.g*self.D**3)*self.rho_d/(self.Ga*self.mu_r);    
        self.mu_f = self.mu_d*self.mu_r
        self.rho_f = self.rho_d*self.rho_r
        self.sig = abs(1-self.rho_r)*self.rho_d*self.g*self.D**2/self.Bo
        # Geometrical parameters
        self.DP = self.Ls /self.N           # size of a unit cell
        self.RND = 0.49*(self.DP-self.D)    # maximum displacement of the buuble from the center of the unit cell
        # numerical parameters
        self.LEVEL =int( math.ceil(np.log(self.ndcmin*self.Ls/self.D)/np.log(2.))  )             # number of cells per diameters of a bubble
        self.ndc = 2**self.LEVEL/self.Ls*self.D
        
    def run(self):
        os.system('make clean')
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir) 
        print('\n')
        print('Study on ',self.para_name,self.parameters_list)
        print('in ',self.namedir)
        print('\n')
        for P in self.parameters_list:
            os.system('rm -rf '+self.name_of_C_file)          
            self.set_parameters_dictionnary(P)
            self.print_parameters()
            print(os.listdir())
            if self.nProc > 1:
                os.system('CC99="mpicc -D_MPI='+str(int(self.nProc)) +'"'+ ' make '+self.name_of_C_file+'.tst')
            else:    
                os.system('make '+self.name_of_C_file+'.tst')
            self.print_parameters_after_sim()
            self.cp_dir_and_files(P)
    
    def set_parameters_dictionnary(self,parameter_value):
        self.update_parameters_dictionnary()
        self.parameters[self.para_name] = parameter_value
        self.calcul_parameters()
        self.update_parameters_dictionnary()
        self.parameters[self.para_name] = parameter_value
    
    def calcul_parameters(self):
        if self.parameters['dimension'] == 2:                                                 
            self.Nb = self.parameters['N']**2
            self.Ls = self.parameters['N']*self.parameters['D']*math.sqrt(math.pi/(4.*self.parameters['PHI']))   
        if self.parameters['dimension'] == 3:                                                  
            self.Nb = self.parameters['N']**3               
            self.Ls = self.parameters['N']*self.parameters['D']*math.pow(math.pi/(6.*self.parameters['PHI']),1./3.)   # size of the domain
        self.mu_d = math.sqrt(abs(1-self.parameters['rho_r'])*self.parameters['rho_r']*self.parameters['g']*self.parameters['D']**3)*self.parameters['rho_d']/(self.parameters['Ga']*self.parameters['mu_r']);    
        self.mu_f = self.mu_d*self.parameters['mu_r']
        self.rho_f = self.parameters['rho_d']*self.parameters['rho_r']
        self.sig = abs(1-self.parameters['rho_r'])*self.parameters['rho_d']*self.parameters['g']*self.parameters['D']**2/self.parameters['Bo']
        # parameters['G']eometrical parameters
        self.DP = self.Ls /self.parameters['N']           # size of a unit cell
        self.RND = 0.49*(self.DP-self.parameters['D'])    # maximum displacement of the buuble from the center of the unit cell
        # numerical parameters
        self.LEVEL =int( math.ceil(np.log(self.parameters['ndcmin']*self.Ls/self.parameters['D'])/np.log(2.)))             # number of cells per diameters of a bubble
        self.ndc = 2**self.LEVEL/self.Ls*self.parameters['D']
        self.SAVESTEP = self.parameters['TMAX']
        
    def update_parameters_dictionnary(self):
        self.parameters = {
            "para_name"         :self.para_name,
            "namestud"          :self.namestud,
            "name_of_C_file"    :self.name_of_C_file,
            "MPI"               :self.MPI,
            "nProc"             :self.nProc,
            "Bo"                :self.Bo,              # Bond number
            "Ga"                :self.Ga,                # Galileo number
            "rho_r"             :self.rho_r,      # ratio of rho
            "mu_r"              :self.mu_r,       # ratio between mu
            "PHI"               :self.PHI,             # concentration phi
            "Nb"                :self.Nb,               # number of bubbles (Must be a perfect square or it will be rounded) 
            "ndc"               :self.ndc,               # number of cells per diameters of a bubble
            "dtprint"           :self.dtprint,          # writing time step
            "dimension"         :self.dimension,      
            "g"                 :self.g,
            "D"                 :self.D,
            "rho_d"             :self.rho_d,
            "N"                 :self.N,
            "Ls"                :self.Ls,
            "mu_d"              :self.mu_d,
            "mu_f"              :self.mu_f,
            "LEVEL"             :self.LEVEL,
            "rho_f"             :self.rho_f,
            "sig"               :self.sig,
            "DP"                :self.DP,
            "r1"                :self.r1,
            "RND"               :self.RND,
            "TMAX"              :self.TMAX,
            "FIRSTSTEP"         :self.FIRSTSTEP,
            "nopara"            :self.nopara,
            "SAVESTEP"          :self.SAVESTEP,
            "MOVIES"            :self.MOVIES,
            "ndcmin"            :self.ndcmin
        }
        
    def print_parameters(self):
        with open('parameters.h','w') as the_file:
            for key,value in self.parameters.items():
                if(type(value) != str and type(value) != bool  ):
                    the_file.write('#define '+key+' '+str(round(value,5))+'\n')
                    
    def print_parameters_after_sim(self):
        with open(self.name_of_C_file+'/parameters.h','w') as the_file:
            for key,value in self.parameters.items():
                if(type(value) != str and type(value) != bool  ):
                    the_file.write('#define '+key+' '+str(round(value,5))+'\n')
        with open(self.name_of_C_file+'/parameters.py','w') as the_file:
            the_file.write('parameters = {\n')
            for key,value in self.parameters.items():
                if(type(value) != str):
                    the_file.write('"'+key+'":'+str(value)+',\n')
                else:    
                    the_file.write('"'+key+'":"'+str(value)+'",\n')
            the_file.write('}')
        
                    
    def cp_dir_and_files(self,P):
        namedir_simulation = self.namedir+self.para_name+'_'+str(P)
        try:
            shutil.rmtree(namedir_simulation)
            os.mkdir(namedir_simulation)
        except:
            os.mkdir(namedir_simulation)
        os.system('cp -r --backup=simple '+self.name_of_C_file+'/*'+' '+namedir_simulation)
        print(os.listdir(''+self.name_of_C_file+'/.'))
        