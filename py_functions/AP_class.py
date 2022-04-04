import os
MYLIBS = os.getenv('MYLIBS')
import re 
import shutil
import numpy as np
import math
import sys
import pandas as pd
from py_functions.generation import *
#sync
class Analyse_Parametrique():
    """Permet de gerer les parametres durant une etude parametrique"""
    def __init__(self,para_name,parameters_list,namestud):
        """This is the initialisation of the AP class it.
        We set up all the defaults parameters that we might need.

        Args:
            para_name (str): This is ne mane of the parameters that will take the values of the parameters_list
            parameters_list (list): values of the parameters that will vary along the study
            namestud (str): directory where the study will be stored (the directory will be located in the self.resultsdir directory)
        """
        ######### Parameters containing the PATH of all usefull directories ##################
        self.resultsdir = "../results/"
        self.bashDir = MYLIBS + "/bash_functions/"
        self.namedir 	= self.resultsdir+namestud.strip()+'/'
        self.para_name 	=  para_name
        self.parameters_list = parameters_list
        self.parameters = {}
        self.ghostparameter = 0
        ##########################################################################:
        ############## Parameters linked to the simulation on OpenFoam ###########:
        ##########################################################################:
        ######## Fixed parameters (meaning that all the others depend on those) ##:
        ############# Dimmensionneless ##################
        self.Re 	= 25
        self.ksi 	= 3
        self.Theta 	= 45
        self.ThetaU     = 0
        self.r 		= 0.5
        self.U 		= 1
        self.Uf 	= 1
        self.rho 	= 1000
        ############# Numerical parameters ####################################
        self.ndc 	= 30 # number of cells per diameter of cylinders
        self.Raffinement_de_surface 		    = 2   # LEVEL of refinement at the surface of the cyl
        self.GapLevel                           = 2   #LEVEL of refinement max to close a gap
        self.Raffinement_de_region 		        = self.Raffinement_de_surface #refinement around a region
        self.ER     = self.Raffinement_de_surface +1 # LEVEL of refinement of the edges
        self.add_layers 			            = 0  #1 if it s true 0 for false
        self.wakezone 				            = 1     # refin the wake zone of true 
        self.nCellsBetweenLevels	            = 2     
        self.xblockmeshrelatif                  = 2  # ratio between X/Y length of the domain
        self.ratio_ZY                           = 1  # ratio between Z/Y length of the domain
        ## Parametre de fvSolution # 
        self.Relaxcoef 	= 0.7
        self.AutoRelaxcoef = 0
        self.nNonOrthogonalCorrectors = 2
        self.EndTime = 2000
        self.ptol = 1e-6
        self.Utol = 1e-7
        self.pres = 1e-5
        self.Ures = 1e-6
        ########## Parametres qui sont fonction des parametres ci dessus ########
        #info sur le cylindre
        self.length    	=self.ksi*self.r*2
        #box lenght
        self.HowManyD   = 20
        self.HowManyL   = 5
        self.x         	=max(self.HowManyD*self.r*2, self.length*self.HowManyL)
        self.x2         =self.x*self.xblockmeshrelatif
        self.y         	=max(self.HowManyD*self.r*2, self.length*self.HowManyL)
        self.z         	=self.y  * self.ratio_ZY
        
        ############ for the cylinder #############
        
        self.xCenter   	=self.x/2
        self.yCenter   	=self.y/2
        self.zCenter   	=self.z/2
        
        #Parametre materiaux 
        #couche limite 
        self.MeshForHightRe = 0
        self.Boundary_layer_lenght 		= self.r*2/(5*math.sqrt(self.Re))
        self.Howmany_cells_per_layers = 5
        self.ndcpd_couche_limite=self.Howmany_cells_per_layers*math.sqrt(self.Re)
        self.Maillage_de_fond = int(self.x2*self.ndc / (2**self.Raffinement_de_surface*2*self.r))
        self.Maillage_de_fondy =int(self.Maillage_de_fond * 1./self.xblockmeshrelatif)
        self.Maillage_de_fondz =int(self.Maillage_de_fond * 1./self.xblockmeshrelatif * self.ratio_ZY)
        self.ReMesh 	        = 1
        self.Auto_refinement    = [0,100]
        self.whichMesh          = 1
        
        ######### For continuing study  ######
        self.REV = Six_pack()
        self.REV.noc = 10
        self.REV.rmean = self.r
        self.REV.lengthmean = self.length
        self.REV.sizex = self.x2
        self.REV.sizey = self.y
        self.REV.sizez = self.z
        
        self.doREV = 1
        self.mPatchs = 0
        self.DeltaP = 1
        #############Parametre de la BC modifie #######################
        self.Relax = 1
        self.UtolBC = 0.01
        
        ########### write only last step of forces ###################
        self.wrinteAtEnd = 0
        self.wrinteEachTimeStep = self.EndTime
        self.TimePimple = 10
        self.printStepPimple = 0.1
        ########### For rotating cylinders ###################
                  
        if self.U != 0:
            self.nu 	= self.r*2*self.U/self.Re
        else:
            self.nu = self.r*2*self.Uz/self.Re
        self.Re_omega = 1.
        self.omega = self.Re_omega*self.nu/(self.length*self.r)
        self.Uz = 0 
        self.Re_omega_para = 0.
        self.omega_para = self.Re_omega_para*self.nu/(self.r*self.r*2)
        
        self.eta       	=self.nu*self.rho

    def run(self):
        """This function launch the parametrical simulation one by one.
        """
        os.system(self.bashDir+'Allclean')
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir) 
        for P  in self.parameters_list:
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            os.system(self.bashDir+'Allclean')
            self.eMeshGenerator()
            self.print_parameters_for_OF()
            os.system(self.bashDir+'Allrun '+str(self.whichMesh))
            self.cp_dir_and_files(P)
            self.parameters_for_Studies(P)

    def run_first_step(self):
        """This function run only the first step of a study.
        """
        os.system(self.bashDir+'Allclean')
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        # run the first step
        P  =self.parameters_list[0]
        self.set_parameters_dictionnary(P)
        print(self.para_name+' = '+str(self.parameters[self.para_name]))
        self.eMeshGenerator()
        self.print_parameters_for_OF()
        os.system(self.bashDir+'Allrun '+str(self.whichMesh))
        self.cp_dir_and_files(P)
        self.parameters_for_Studies(P)
        # run parametric

    def run_parametric(self):
        """This function take the step n-1 to start the n-th simulation."""
        a = self.parameters_list[1:]
        b = self.parameters_list[:-1]
        if len(self.parameters_list) == 1:
            a  = (self.parameters_list*2)[1:]
            b  = (self.parameters_list*2)[:-1]
        for P,Pmoins1  in zip(a,b):
            os.system(self.bashDir+'Allclean')
            os.system('rm -f log.*')
            self.cp_first_step_here(Pmoins1)
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.print_parameters_for_OF()
            os.system(self.bashDir+'AllrunNotrestore0Dir  '+str(self.whichMesh))
            self.cp_dir_and_files(P)
            self.parameters_for_Studies(P)

            
    def run_pimple_parametric(self):
        """Same as run_parametric but with unsteady the solver pimpleFoam"""
        a = self.parameters_list[1:]
        b = self.parameters_list[:-1]
        if len(self.parameters_list) == 1:
            a  = (self.parameters_list*2)[1:]
            b  = (self.parameters_list*2)[:-1]
            
        for P,Pmoins1  in zip(a,b):
            os.system(self.bashDir+'Allclean')
            os.system('rm -f log.*')
            self.cp_first_step_here(Pmoins1)
            self.wrinteEachTimeStep = self.EndTime
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.print_parameters_for_OF()
            self.printPatch()
            os.system(self.bashDir+'Allrun2bis')
            dirs = pd.Series(os.listdir('processor0')) 
            self.EndTime  = float(dirs[dirs.str.isdigit()].max()) + self.TimePimple
            self.wrinteEachTimeStep = self.printStepPimple
            self.UtolBC = 1
            self.Relaxcoef = 1
            self.set_parameters_dictionnary(P)
            self.print_parameters_for_OF()
            os.system(self.bashDir+'Allrun3')
            self.cp_dir_and_files(P)
            self.parameters_for_Studies(P)

######################### The point wise runs functions ################
    def run_PW(self):
        """Run a simulation with a pre created mesh on pPointWise for example"""
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)

        for P  in self.parameters_list:
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.print_parameters_for_OF()
            os.system('rm -f log.*')
            os.system(self.bashDir+'AllrunExceptTheMesh')
            self.mv_dir_and_files(P)
            self.parameters_for_Studies(P)
    
    def run_PW_first_step(self):
        """Run the first step of a study where the mesh is already created
        """
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        P  =self.parameters_list[0]
        self.set_parameters_dictionnary(P)
        print(self.para_name+' = '+str(self.parameters[self.para_name]))
        self.print_parameters_for_OF()
        os.system('rm -f log.*')
        os.system(self.bashDir+'AllrunExceptTheMesh')
        self.cp_dir_and_files(P)
        self.parameters_for_Studies(P)


################# REV runs functions ##############################
    def run_REV(self):
        """Launch the parametrical study for REV"""
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        for P in self.parameters_list:
            if self.REV.DoGen:
                self.REV.Cyls = []
            os.system(self.bashDir+'Allclean')
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.eMeshRVE()
            self.Uf = self.U/(1.-self.REV.phi)
            print("mean fluid vel :",self.Uf)
            self.set_parameters_dictionnary(P)
            self.print_parameters_for_OF()
            os.system(self.bashDir+'Allrun1 '+str(self.whichMesh))
            self.printPatch()
            self.forcesFunc('p','U')
            os.system(self.bashDir+'Allrun1bis '+str(self.whichMesh))
            os.system(self.bashDir+'Allrun2 '+str(self.whichMesh))
            self.Cyls_for_Studies(P)
            self.cp_dir_and_files(P)
            self.parameters_for_Studies(P)
            
    def run_REV_pimple(self):
        """Launch the parametrical study for REV with pimpleFoam"""
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        for P in self.parameters_list:
            if self.REV.DoGen:
                self.REV.Cyls = []
            os.system(self.bashDir+'Allclean')
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.eMeshRVE()
            self.Uf = self.U/(1.-self.REV.phi)
            print("mean fluid vel :",self.Uf)
            self.wrinteEachTimeStep = self.EndTime
            self.set_parameters_dictionnary(P)
            self.print_parameters_for_OF()
            os.system(self.bashDir+'Allrun1 '+str(self.whichMesh))
            self.printPatch()
            self.forcesFunc('p','U')
            os.system(self.bashDir+'Allrun1bis '+str(self.whichMesh))
            os.system(self.bashDir+'Allrun2bis')
            dirs = pd.Series(os.listdir('processor0')) 
            self.EndTime  = float(dirs[dirs.str.isdigit()].max()) + self.TimePimple
            self.wrinteEachTimeStep = self.printStepPimple
            self.UtolBC = 1
            self.Relaxcoef = 1
            self.set_parameters_dictionnary(P)
            self.print_parameters_for_OF()
            os.system(self.bashDir+'Allrun3')
            self.Cyls_for_Studies(P)
            self.cp_dir_and_files(P)
            self.parameters_for_Studies(P)
            
    def run_Mesh(self):
        """Generate a mesh for REV study. 
        The aim of this function is to generate some meshes that works fine.
        Then in a second study one can run simulations with thoses mesh using run_parametric()
        """
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        P = self.parameters_list[0]
        if self.REV.DoGen:
            self.REV.Cyls = []
        os.system(self.bashDir+'Allclean')
        self.set_parameters_dictionnary(P)
        print(self.para_name+' = '+str(self.parameters[self.para_name]))
        self.eMeshRVE()
        self.Uf = self.U/(1.-self.REV.phi)
        print("mean fluid vel :",self.Uf)
        self.wrinteEachTimeStep = self.EndTime
        self.EndTime = 2
        self.set_parameters_dictionnary(P)
        self.print_parameters_for_OF()
        os.system(self.bashDir+'Allrun1')
        self.printPatch()
        self.forcesFunc('p','U')
        os.system(self.bashDir+'Allrun1bis')
        os.system(self.bashDir+'Allrun2bis')
        self.Cyls_for_Studies(P)
        self.cp_dir_and_files(P)
        self.parameters_for_Studies(P)
            

    def run_REV_first_step(self):
        """Run the first step of a REV study.
        It is adviced to use run_mesh, since the meshes often don't work.
        """
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        os.system(self.bashDir+'Allclean')
        P  =self.parameters_list[0]
        self.set_parameters_dictionnary(P)
        print(self.para_name+' = '+str(self.parameters[self.para_name]))
        self.eMeshRVE()
        self.Uf = self.U/(1.-self.REV.phi)
        print("mean fluid vel :",self.Uf)
        self.set_parameters_dictionnary(P)
        self.print_parameters_for_OF()
        os.system(self.bashDir+'Allrun1 '+str(self.whichMesh))
        self.printPatch()
        self.forcesFunc('p','U')
        os.system(self.bashDir+'Allrun1bis '+str(self.whichMesh))
        os.system(self.bashDir+'Allrun2 '+str(self.whichMesh))
        self.Cyls_for_Studies(P)
        self.cp_dir_and_files(P)
        self.parameters_for_Studies(P)

###########################################################################

    def getStudyParameters(self): 
        sys.path.append('./'+self.namedir)
        import parameters_for_this_study as pstud
        ###### Update the parameters liked to the mesh ########
        self.ksi = pstud.parameters['ksi']
        self.HowManyD = pstud.parameters['HowManyD']
        self.HowManyL = pstud.parameters['HowManyL']
        self.ratio_ZY = pstud.parameters['ratio_ZY']
        self.r = pstud.parameters['r']
        try:
            self.ndc = pstud.parameters['ndc']
        except:
            self.ndc = pstud.parameters['nombre_de_cellules_par_diametre']
        try:
            self.whichMesh = pstud.parameters['whichMesh']
        except:
            1

        self.Raffinement_de_surface = pstud.parameters['Raffinement_de_surface']

    def set_parameters_dictionnary(self,parameter_value):
        self.update_parameters_dictionnary()
        self.parameters[self.para_name] = parameter_value
        self.calcul_parameters(self.parameters)
        self.update_parameters_dictionnary()
        self.parameters[self.para_name] = parameter_value

    def calcul_parameters(self,parameters):
    
        ########## Parametres fonction des parametres ci dessus ######## 
        ############## if self.para_name == 'Re' use this formula : ect############
        #info sur le cylindre
        self.length    =parameters['ksi']*parameters['r']*2
        #box lenght
        self.x         =max(parameters['HowManyD']*parameters['r']*2,self.length*parameters['HowManyL'])
        self.x2         =self.x*self.xblockmeshrelatif
        #Parametre materiaux 
        if self.whichMesh == 4:
            self.U = 0.
            self.Uz = 1.
        else:
            self.Uz = 0.
            
        if self.parameters['U'] != 0:
            self.nu = parameters['r']*2*parameters['U']/parameters['Re']
        else:
            self.nu = parameters['r']*2*parameters['Uz']/parameters['Re']
        #couche limite 
        self.Boundary_layer_lenght = parameters['r']*2/(5*math.sqrt(parameters['Re']))
        self.ndcpd_couche_limite=self.Howmany_cells_per_layers*math.sqrt(parameters['Re'])
        
        if self.MeshForHightRe == 1:
            self.ndc = self.ndcpd_couche_limite
        #Maillage a la couche limite 
        self.Maillage_de_fond = int(self.x2*self.parameters['ndc'] / (2**self.parameters['Raffinement_de_surface']*2*self.parameters['r']))

        self.Raffinement_de_region = self.parameters['Raffinement_de_surface'] 
        #self.ER = self.parameters['Raffinement_de_surface'] + 1

        self.Maillage_de_fondy =int(self.Maillage_de_fond * 1./self.xblockmeshrelatif)
        
        self.Maillage_de_fondz =int(self.Maillage_de_fond * 1./self.xblockmeshrelatif  * self.parameters['ratio_ZY'])
        ############# Parametres qui dependent des para qui viennent detre calcule########
        self.y         =max(self.x, self.length*parameters['HowManyL'])
        self.z         =self.y  * self.parameters['ratio_ZY']
        self.xCenter   =self.x/2
        self.yCenter   =self.y/2
        self.zCenter   =self.z/2
        #Parametre materiaux 
        self.eta       =self.nu*parameters['Rho']
        if self.AutoRelaxcoef == 1:
            if self.parameters['Re'] < 5: 
                self.Relaxcoef = 0.95
            if self.parameters['Re'] <= 0.03:
                self.Relaxcoef = 0.98
            if self.parameters['Re'] > 1 and self.parameters['Re'] <= 10: 
                self.Relaxcoef = 0.85
            if self.parameters['Re'] > 10 and self.parameters['Re'] <= 80: 
                self.Relaxcoef = 0.65
            if self.parameters['Re'] > 80 : 
                self.Relaxcoef = 0.65
        ########### Parametes du REV a re calculer ######################
        self.REV.noc = self.parameters["noc"]
        self.REV.rmean = self.parameters['r']
        self.REV.lengthmean = self.length
        self.REV.sizex = self.x2
        self.REV.sizey = self.y
        self.REV.sizez = self.z
        self.forces = []
        if self.whichMesh == 3 or self.whichMesh == 4:
            self.xCenter = 0
            self.yCenter = 0
            
        ########### For rotating cylinders ###################
        self.Re_omega = self.parameters['Re_omega']
        self.Re_omega_para = self.parameters['Re_omega_para']
        self.omega = self.Re_omega*self.nu/(self.length*self.parameters['r'])
        self.omega_para = self.Re_omega_para*self.nu/(self.parameters['r']*self.parameters['r']*2)
       
            

    
    def update_parameters_dictionnary(self):
        self.parameters = {
            "Uz"        :self.Uz,
            "omega"     :self.omega,
            "Re_omega"  :self.Re_omega,
            "Re_omega_para"  :self.Re_omega_para,
            "omega_para"     :self.omega_para,
            "Uf"        :self.Uf,
            "UtolBC"    :self.UtolBC,
            "Relax"     :self.Relax,
            "eta"       :self.eta,
            "GapLevel"  :self.GapLevel,
            "whichMesh" :self.whichMesh,
            "Theta"     :self.Theta,
            "ThetaU"    :self.ThetaU,
            "length"    :self.length,
            "r"         :self.r,
            #box len,
            "HowManyD"  :self.HowManyD,
            "HowManyL"  :self.HowManyL,
            "x"         :self.x,
            "x2"        :self.x2,
            "y"         :self.y ,
            "z"         :self.z,
            "xCenter"   :self.xCenter,
            "yCenter"   :self.yCenter,
            "zCenter"   :self.zCenter,
            #level of eMesh rafineme,
            "ER"        :self.ER,
            "Re"        :self.Re,
            "nu"        :self.nu,
            "Rho"       :self.rho,
            "U"         :self.U,
            "ksi"       :self.ksi,
            "TimePimple":self.TimePimple,
            "Maillage_de_fond":self.Maillage_de_fond,
            "Boundary_layer_lenght":self.Boundary_layer_lenght,
            "Raffinement_de_surface":self.Raffinement_de_surface,
            "ndc":self.ndc,
            "Raffinement_de_region":self.Raffinement_de_region,
            "para_name":self.para_name,
            "namedir":self.namedir,
            "add_layers":self.add_layers,
            "Relaxcoef":self.Relaxcoef,
            "nNonOrthogonalCorrectors":self.nNonOrthogonalCorrectors,
            "EndTime":self.EndTime,
            "ghostparameter":self.ghostparameter,
            "wakezone":self.wakezone,
            "nCellsBetweenLevels":self.nCellsBetweenLevels,
            "Maillage_de_fondy":self.Maillage_de_fondy,
            "Maillage_de_fondz":self.Maillage_de_fondz,
            "ptol":self.ptol,
            "Utol":self.Utol,
            "pres":self.pres,
            "Ures":self.Ures,
            "MeshForHightRe":self.MeshForHightRe,
            "ratio_ZY":self.ratio_ZY,
            "Mx":self.REV.Matiere[0],
            "My":self.REV.Matiere[1],
            "Mz":self.REV.Matiere[2],
            "T11":self.REV.T[0],
            "T22":self.REV.T[1],
            "T33":self.REV.T[2],
            "pRefx":self.REV.pRef[0],
            "pRefy":self.REV.pRef[1],
            "pRefz":self.REV.pRef[2],
            "noc":self.REV.noc,
            "DeltaP"    :self.DeltaP,
            "Shape"     :self.REV.shape,
            "phi"       :self.REV.phi,
            "wrinteEachTimeStep":self.wrinteEachTimeStep
        }

    def eMeshGenerator(self):
        """Create all the vectors of the eMesh"""
        t =math.radians(self.parameters['Theta'])
        Ocyl = [self.parameters['xCenter'],self.parameters['yCenter'],self.parameters['zCenter']]
        # Coordonne d'un point M du cercle du haut parametre avec l'angle gamma 
        gamma = list(np.linspace(0,math.pi*2,5000))


        for numero_du_cercle,l in [0,self.parameters['length']],[1,-self.parameters['length']] :
            M_x = [Ocyl[0]+l/2*math.cos(t)-self.parameters['r']*math.cos(g)*math.sin(t)for g in gamma]
            M_y = [Ocyl[1]+l/2*math.sin(t)+self.parameters['r']*math.cos(g)*math.cos(t) for g in gamma]
            M_z = [Ocyl[2] + self.parameters['r']*math.sin(g) for g in gamma]
            M = list(zip(M_x,M_y,M_z))
            name="cercle"+str(numero_du_cercle)
            self.eMesh(M,name,1)
        self.eMeshForSnappy()
    
    def eMeshRVE(self):
        if self.doREV == 1:
            self.REV.run()
                
        os.system('rm -rf constant/triSurface/*')
        for Cyl,i in zip(self.REV.Cyls,range(len(self.REV.Cyls))):
            for key in Cyl[0].keys:
                if Cyl[0].M[key] != []:
                    name = 'cercle'+'_'+str(Cyl[1])+'_'+key
                    self.eMesh(Cyl[0].M[key],name,0)

        ##### for sphere #####
        for Sphere,i in zip(self.REV.Spheres,range(len(self.REV.Spheres))):
            for key in Sphere[0].keys:
                if Sphere[0].M[key] != [] :
                    name = 'cercle'+'_'+str(Sphere[1])+'_'+key
                    self.eMesh(Sphere[0].M[key],name,0)

        self.CylsForSnappy()
        self.eMeshForSnappy()
    
    def eMesh(self,M,name,loop):
        """print the eMesh files from a set of points M"""
        eMesh = []    
        eMesh.append("/*--------------------------------*- C++ -*----------------------------------*\\")
        eMesh.append("| =========                 |                                                  |")
        eMesh.append("| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |")
        eMesh.append("|  \\\    /   O peration     | Version:  plus                                  |")
        eMesh.append("|   \\\  /    A nd           | Website:  www.openfoam.com                      |")
        eMesh.append("|    \\\/     M anipulation  |                                                 |")
        eMesh.append("\*---------------------------------------------------------------------------*/")
        eMesh.append("FoamFile")
        eMesh.append("{")
        eMesh.append("\t"+"version"+"\t"+"\t"+"2.0;")
        eMesh.append("\t"+"format"+"\t"+"\t"+"ascii;")
        eMesh.append("\t"+"class"+"\t"+"\t"+"featureEdgeMesh;")
        eMesh.append("\t"+"location"+"\t"+'"'+'constant/triSurface'+'";')
        eMesh.append("\t"+"object"+"\t"+"\t"+name+'.eMesh;')
        # eMesh.append("\t"+"object"+"\t"+"\t"+"cyl"+repr(N)+".eMesh;")
        eMesh.append("}")
        eMesh.append("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"+"\n"+"\n")    
        eMesh.append("// ** This file was created by the python code AP_class.py")
        # plot des points
        eMesh.append(len(M))
        eMesh.append('(')
        for x,y,z in M:
            eMesh.append('('+repr(x)+'\t'+repr(y)+'\t'+repr(z)+')')
        eMesh.append(")"+"\n"+"\n")  
        #plot des edgs
        if loop:
            eMesh.append(len(M))
        else :
            eMesh.append(len(M)-1)
        eMesh.append('(')
        for i in range(0,len(M)-1): #1st circle
            eMesh.append('('+repr(i)+'\t'+repr(i+1)+')')
        if loop:
            eMesh.append('('+repr(len(M)-1)+'\t'+repr(0)+')')
        eMesh.append(')')
        eMesh.append("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"+"\n")
        eMesh.append("// ************************************************************************* //")
        np.savetxt('constant/triSurface/'+name+'.eMesh',eMesh, fmt='%s',delimiter=" ")

    def eMeshForSnappy(self):
        alleMeshnames = filter(lambda x : x[-5:]=="eMesh" ,   os.listdir('constant/triSurface'))
        print(alleMeshnames  )
        eMeshforSnappy = []
        for eMeshname in alleMeshnames:
            # import parameters
            eMeshforSnappy.append("{")
            eMeshforSnappy.append('\t file '+'"'+eMeshname+'";')
            # eMeshforSnappy.append("\t levels (("+str(self.r)+"    "+str(self.parameters['ER'])+"));")
            eMeshforSnappy.append("\t level   "+str(self.parameters['ER'])+";")
            eMeshforSnappy.append("}\n") 
        np.savetxt('system/eMeshforSnappy',eMeshforSnappy, fmt='%s',delimiter=" ")

    def CylsForSnappy(self):
        import numpy as np
        # import parameters
        cylforsnap = []
        for Cyl,i in zip(self.REV.Cyls,range(len(self.REV.Cyls))):
            name = Cyl[1]
            Cyl = Cyl[0]
            cylforsnap.append("air_to_sphere"+str(name))
            cylforsnap.append("{")
            cylforsnap.append('\ttype searchableCylinder;')
            cylforsnap.append("\tpoint1 ("+ str(Cyl.p1[0])+" "+str(Cyl.p1[1])+" "+str(Cyl.p1[2])+");")
            cylforsnap.append("\tpoint2 ("+str(Cyl.p2[0])+" "+str(Cyl.p2[1])+" "+str(Cyl.p2[2])+ ");")
            cylforsnap.append ("\tradius "+ str(Cyl.r)+";")
            cylforsnap.append("\tname air_to_sphere"+str(name)+";")
            cylforsnap.append("}\n") 
 
        for Sphere,i in zip(self.REV.Spheres,range(len(self.REV.Spheres))):
            name = Sphere[1]
            Sphere = Sphere[0]
            cylforsnap.append("air_to_sphere"+str(name))
            cylforsnap.append("{")
            cylforsnap.append('\ttype searchableSphere;')
            cylforsnap.append("\tcentre ("+ str(Sphere.OP[0])+" "+str(Sphere.OP[1])+" "+str(Sphere.OP[2])+");")
            cylforsnap.append ("\tradius "+ str(Sphere.r)+";")
            cylforsnap.append("\tname air_to_sphere"+str(name)+";")
            cylforsnap.append("}\n") 
        np.savetxt('system/cylforsnap',cylforsnap, fmt='%s',delimiter=" ")

    def printPatch(self):
        # import parameters
        patch = []
        nameP1 = open('patch_list.txt','r').readlines()
        if self.mPatchs == 1:
            for i in range(self.REV.noc):
                names = [name.rstrip("\n") for name in nameP1 if name.split('_')[2].rstrip("\n") =='sphere'+str(i)]
                name = ''
                for n in names:
                    name = name +' '+n
                patch.append('{')
                patch.append('\tname Cylinder_'+str(i)+';')
                patch.append('\tpatchInfo')
                patch.append('\t{')
                patch.append('\t  \t type wall;')
                patch.append('\t}')
                patch.append('\tconstructFrom patches;')
                patch.append('\tpatches ('+name+');')
                # patch.append('\tpatches (air_to_sphere'+str(i)+'.*);')
                patch.append('}')                          
        np.savetxt('system/patch',patch, fmt='%s',delimiter=" ")

    def forcesFunc(self,p,U):
        # import parameters
        nameP1 = open('patch_list.txt','r').readlines()
        names = [name.rstrip("\n") for name in nameP1]
        #creating the list for the true Cyl no yet merged
        Ids  = [name[17:] for name in names]
        self.REV.VeryRealCyl = [Cyl for Cyl in self.REV.Cyls if Cyl[1] in Ids]
        self.REV.VeryRealSphere = [Sphere for Sphere in self.REV.Spheres if Sphere[1] in Ids]
        if self.mPatchs == 0:
            for Cyl,i,name in zip(self.REV.VeryRealCyl,range(len(self.REV.VeryRealCyl)),names):
                self.forces.append('forces_'+Cyl[1])
                self.forces.append('{')
                self.forces.append('\ttype            FirstMoment;')
                self.forces.append('\tlibs            (libFirstMoment);')
                self.forces.append('\twriteControl    timeStep;')
                self.forces.append('\twriteInterval    $wrinteEachTimeStep;')
                self.forces.append('\texecuteControl    timeStep;')
                self.forces.append('\texecuteInterval   $wrinteEachTimeStep;')
                self.forces.append('\tpatches         ('+name+');')
                self.forces.append('\tp              '+p+';')                        
                self.forces.append('\tU              '+U+';')                        
                self.forces.append('\tlog             true;')                        
                self.forces.append('\trho             rhoInf;')  
                self.forces.append('\tCofR            ('+str(Cyl[0].OP[0])+" "+str(Cyl[0].OP[1])+" "+str(Cyl[0].OP[2])+');')                        
                self.forces.append('\trhoInf          $Rho;')       
                if Cyl[0].pts['Ixb']['X'] != []:
                    self.forces.append('\tpRef           $jump;')             
                self.forces.append('}')                      
            for Sphere,i,name in zip(self.REV.VeryRealSphere,range(len(self.REV.VeryRealSphere)),names):
                self.forces.append('forces_'+Sphere[1])
                self.forces.append('{')
                self.forces.append('\ttype            FirstMoment;')
                self.forces.append('\tlibs            (libFirstMoment);')
                self.forces.append('\twriteControl    timeStep;')
                self.forces.append('\twriteInterval    $wrinteEachTimeStep;')
                self.forces.append('\texecuteControl    timeStep;')
                self.forces.append('\texecuteInterval   $wrinteEachTimeStep;')
                self.forces.append('\tpatches         ('+name+');')
                self.forces.append('\tp              '+p+';')                        
                self.forces.append('\tU              '+U+';')                        
                self.forces.append('\tlog             true;')                        
                self.forces.append('\trho             rhoInf;')  
                self.forces.append('\tCofR            ('+str(Sphere[0].OP[0])+" "+str(Sphere[0].OP[1])+" "+str(Sphere[0].OP[2])+');')                        
                self.forces.append('\trhoInf          $Rho;')       
                if Sphere[0].pts['Ixb']['X'] != []:
                    self.forces.append('\tpRef           $jump;')             
                self.forces.append('}')                      
          
        np.savetxt('system/forces',self.forces, fmt='%s',delimiter=" ")

    def print_parameters_for_OF(self):
        """print the parameter C++ file for OpenFoam"""
        parametersOF = []  
        parametersOF.append("/*--------------------------------*- C++ -*----------------------------------*\\")
        parametersOF.append("| =========                 |                                                 |")
        parametersOF.append("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |")
        parametersOF.append("|  \\    /   O peration     | Version:  plus                                  |")
        parametersOF.append("|   \\  /    A nd           | Website:  www.openfoam.com                      |")
        parametersOF.append("|    \\/     M anipulation  |                                                 |")
        parametersOF.append("\*---------------------------------------------------------------------------*/")
        parametersOF.append("FoamFile")
        parametersOF.append("{")
        parametersOF.append("    version     2.0;")
        parametersOF.append("    format      ascii;")
        parametersOF.append('    arch        "LSB;label=32;scalar=64";')
        parametersOF.append("    class       dictionary;")
        parametersOF.append('    location    "system";')
        parametersOF.append("    object      parametersOF;")
        parametersOF.append("}")
        parametersOF.append("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")
        for key,value in list(self.parameters.items()):# for python2 :: .iteritems() :
            if type(value) != str:
                parametersOF.append(key+"      "+str(value)+";"+"\n")
        np.savetxt('system/parametersOF',parametersOF, fmt='%s',delimiter=" ")

    def cp_dir_and_files(self,P):
        namedir_simulation = self.namedir+str(P)
        try:
            print(os.listdir('postProcessing/.'))
        except:
            print('no postProcess stop all the sim')
            quit()            
        
        try:
            shutil.rmtree(namedir_simulation)
            os.mkdir(namedir_simulation)
        except:
            os.mkdir(namedir_simulation)
        os.system('cp -r --backup=simple postProcessing/*'+' '+namedir_simulation)
        os.system('cp -r --backup=simple log*  '+namedir_simulation)
            # that is for having the name of the last dir
        dirs = os.listdir('.')
        dirs = filter(lambda x : x!='' and x!='0.' , [re.sub("[^0123456789\.]","",dir) for dir in dirs])
        dirs =[dire for dire in dirs if dire in os.listdir('.')]
        dirs.sort(key = lambda x : float(x))#range dans l'ordre croissant 
        lastdir = dirs[-1]	#on prend le last car c'est celui qui nous interresse (enfait on prend tt)
        for dire in dirs:  
                os.system('cp -r '+dire+' '+namedir_simulation+'/')
        os.system('cp -r constant '+namedir_simulation)
        os.system('cp -r system '+namedir_simulation)
        os.system('cp -r Cyls_for_this_study.py '+namedir_simulation)
        # os.system('cp -r 0 '+namedir_simulation)
        
    def mv_dir_and_files(self,P):
        namedir_simulation = self.namedir+str(P)
        try:
            print(os.listdir('postProcessing/.'))
        except FileNotFoundError:
            print('no postProcess stop all the sim')
            quit()  
        try:
            shutil.rmtree(namedir_simulation)
            os.mkdir(namedir_simulation)
        except:
            os.mkdir(namedir_simulation)
        print(os.listdir('postProcessing/.'))
        os.system('mv --backup=simple postProcessing/*'+' '+namedir_simulation)
        os.system('mv --backup=simple log*  '+namedir_simulation)	
            # that is for having the name of the last dir
        dirs = os.listdir('.')
        dirs = filter(lambda x : x!='' and x !='.' and x!='0.' , [re.sub("[^0123456789\.]","",dir) for dir in dirs])
        dirs =[dire for dire in dirs if dire in os.listdir('.')]
        dirs.sort(key = lambda x : float(x))#range dans l'ordre croissant 
        lastdir = dirs[-1]	#on prend le last car c'est celui qui nous interresse (enfait on prend tt)
        for dire in dirs:  
                os.system('mv --backup=simple '+dire+' '+namedir_simulation+'/')
        os.system('cp -r constant '+namedir_simulation)
        os.system('cp -r system '+namedir_simulation)
        self.Cyls_for_Studies(P)
        self.parameters_for_Studies(P)

    def parameters_for_Studies(self,P = ''):
        namedir_simulation = self.namedir+str(P)
        """print the parameter file for OpenFoam post traitement"""
        files = []    
        files.append('parameters = {')
        for key,value in self.parameters.items(): #.iteritems() : for python2
            if type(value) != str:
                files.append('  "'+key+'":'+str(value)+',')
            else:
                files.append('  "'+key+'":'+'"'+str(value)+'"'+',')
        # files.append('  "para_name":'+'"'+str(self.para)+'"'+',')
        # files.append('  "namedir":'+'"'+str(namedir)+'"')
        files.append('}')
        print(namedir_simulation)
        np.savetxt(self.namedir+'/parameters_for_this_study.py',files, fmt='%s',delimiter=" ")
        np.savetxt(namedir_simulation+'/parameters_for_this_study.py',files, fmt='%s',delimiter=" ")
        
    def Cyls_for_Studies(self,P =''):
        """print the Cyl C++ file for OpenFoam"""
        namedir_simulation = self.namedir+str(P)
        files = []    
        files.append('Ids = [')
        for Cyl in self.REV.Cyls :
            files.append('\t'+'"'+  Cyl[1]  +'"'+',')
        files.append(']')
        files.append('pos = [')
        for Cyl in self.REV.Cyls :
            files.append('\t['+str(Cyl[0].OP[0])+","+str(Cyl[0].OP[1])+","+str(Cyl[0].OP[2])+'],')
        files.append(']')
        files.append('ori = [')
        for Cyl in self.REV.Cyls :
            files.append('\t['+str(Cyl[0].e1[0])+","+str(Cyl[0].e1[1])+","+str(Cyl[0].e1[2])+'],')
        files.append(']')
        files.append('dim = [')
        for Cyl in self.REV.Cyls :
            files.append('\t['+str(Cyl[0].r)+","+str(Cyl[0].length)+'],')
        files.append(']')
        # np.savetxt(namedir_simulation+'/Cyls_for_this_study.py',files, fmt='%s',delimiter=" ")
        np.savetxt('Cyls_for_this_study.py',files, fmt='%s',delimiter=" ")

    
    def cp_first_step_here(self,P):
        '''this function take the constant and last step to the case'''
        namedir_simulation = self.namedir+str(P)+'/'
        # getting the last dir 
        dirs = os.listdir(namedir_simulation)
        dirs = pd.Series(os.listdir(namedir_simulation))
        dirs = dirs[dirs.str.isdigit()]
        list(filter(lambda x: 'U' in os.listdir(namedir_simulation+x), dirs.values))
        laststep  = str(max([int(a) for a in dirs[dirs.str.isdigit()]]))
        laststep_file = namedir_simulation + laststep
        polyMesh = namedir_simulation + 'constant/polyMesh/'

        # mv the files here 
        os.system('cp -r '+polyMesh+' constant/.')
        os.system('cp -r '+laststep_file+' '+laststep)
        os.system('rm -rf system/forces')
        os.system('cp '+namedir_simulation+'system/forces system/.')
        os.system('cp '+namedir_simulation+'Cyls_for_this_study.py .')
        # import the parameters related to the mesh
        sys.path.append(namedir_simulation)
        import parameters_for_this_study as pstud
        from importlib import reload
        reload(pstud)
        self.parameters = pstud.parameters
        self.Uf = pstud.parameters['Uf']
        self.REV.Matiere[0] = pstud.parameters["Mx"]
        self.REV.Matiere[1] = pstud.parameters["My"]
        self.REV.Matiere[2] = pstud.parameters["Mz"]
        self.REV.pRef[0] = pstud.parameters["pRefx"]
        self.REV.pRef[1] = pstud.parameters["pRefy"]
        self.REV.pRef[2] = pstud.parameters["pRefz"]
        self.REV.phi = pstud.parameters["phi"]
        self.REV.T[0] = pstud.parameters["T11"]
        self.REV.T[1] = pstud.parameters["T22"]
        self.REV.T[2] = pstud.parameters["T33"]
        self.REV.noc = pstud.parameters["noc"]
        sys.path.remove(namedir_simulation)
    
    def cp_constant_file_here(self,P):        
        '''this function take the constant dir to the case'''
        namedir_simulation = self.namedir+str(P)+'/'
        constant_file = namedir_simulation + 'constant'
        os.system('rm -rf constant/')

        os.system('cp -rf  '+constant_file+' constant')

