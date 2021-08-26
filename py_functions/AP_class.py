import os
import re 
import shutil
import numpy as np
import math
import sys
sys.path.append('../py_functions')
import generation 
#sync
class Analyse_Parametrique():
    """Permet de gerer les parametres durant une etude parametrique"""
    def __init__(self,para_name,parameters_list,namestud):
        #######parametres de l etude parametrique##########
        #le parametre selectionne est aussi celui qui va etre post traite
        self.resultsdir = "../results/"
        self.bashDir = "../myLibs/bash_functions/"
        self.para_name 	=  para_name
        self.namedir 	= self.resultsdir+namestud.strip()+'/'
        # list des valeur du parametre a run
        self.parameters_list = parameters_list
        self.parameters = {}
        self.ghostparameter = 0
        ############ Parametres independant et sans dimension ####################:
        ############# Sans Dim ##################
        self.Re 	= 25
        self.ksi 	= 3
        self.Theta 	= 45
        self.ThetaU     = 0
        ############ Parametres Fixe #############
        self.r 		= 0.1
        self.U 		= 1
        self.Uf 	= 1
        self.rho 	= 1000
        ########### Parametres Independant ##########
        #Parametres du maillage 
        # self.Maillage_de_fond = 70 #par block mesh 
        self.ndc 	= 28
        self.Raffinement_de_surface 		    = 3
        self.GapLevel                           = 3
        self.Raffinement_de_region 		        = self.Raffinement_de_surface
        self.ER     = self.Raffinement_de_surface +1
        self.add_layers 			            = 0 #1 if it s true 0 for false
        self.wakezone 				            = 1
        self.nCellsBetweenLevels	            = 2
        self.xblockmeshrelatif                  = 2 # can't be 0
        self.ratio_ZY                           = 1
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
        self.nu 	= self.r*2*self.U/self.Re
        self.eta       	=self.nu*self.rho
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
        self.REV = generation.Six_pack()
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

    def run(self):
        os.system(self.bashDir+'Allclean')
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir) 
        if self.ReMesh == 1:
            self.run_RM()
        elif self.ReMesh == 0:
            self.run_NM()

    def run_RM(self):
        """this is the run function it run every other function and loop thought it"""
        for P  in self.parameters_list:
            self.set_parameters_dictionnary(P)
            self.refinement_if_needed(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.eMeshGenerator()
            self.print_parameters_for_OF()
            os.system(self.bashDir+'Allclean')
            os.system(self.bashDir+'Allrun '+str(self.whichMesh))
            self.cp_dir_and_files(P)
            self.parameters_for_Studies()
        
    def run_NM(self):
        '''That function won't re do the mesh after all studies'''
        P = self.parameters_list[0]
        self.set_parameters_dictionnary(P)
        self.refinement_if_needed(P)
        print(self.para_name+' = '+str(self.parameters[self.para_name]))
        self.eMeshGenerator()
        self.print_parameters_for_OF()
        os.system(self.bashDir+'Allclean')
        os.system(self.bashDir+'RunMeshOnly ' + str(self.whichMesh))
        for P  in self.parameters_list:
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.print_parameters_for_OF()
            os.system(self.bashDir+'AllrunExceptTheMesh')
            self.mv_dir_and_files(P)
            self.parameters_for_Studies()

    def run_FAOC(self):
        '''That function will run an old cases'''
        os.system(self.bashDir+'Allclean')
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        for P  in self.parameters_list:
            self.cp_first_step_here(P)
            self.cp_constant_file_here(P)
            self.getStudyParameters()
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.print_parameters_for_OF()
            os.system(self.bashDir+'G')
            self.mv_dir_and_files(P)
            self.parameters_for_Studies()
            #then i have to update parameters...

    def run_first_step(self):
        # Run the first step from nothing 
        os.system(self.bashDir+'Allclean')
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        # run the first step
        P  =self.parameters_list[0]
        self.set_parameters_dictionnary(P)
        self.refinement_if_needed(P)
        print(self.para_name+' = '+str(self.parameters[self.para_name]))
        self.eMeshGenerator()
        self.print_parameters_for_OF()
        os.system(self.bashDir+'Allrun '+str(self.whichMesh))
        self.cp_dir_and_files(P)
        self.parameters_for_Studies()
        # run parametric

    def run_parametric(self):
        for P,Pmoins1  in zip(self.parameters_list[1:],self.parameters_list[:-1]):
            self.cp_first_step_here(Pmoins1)
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.print_parameters_for_OF()
            os.system('rm -f log.*')
            os.system(self.bashDir+'AllrunNotrestore0Dir')
            self.mv_dir_and_files(P)
            self.Cyls_for_Studies(P)
            self.parameters_for_Studies()

######################### The point wise runs functions ################
    def run_PW(self):
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)

        for P  in self.parameters_list:
            self.set_parameters_dictionnary(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.print_parameters_for_OF()
            os.system('rm -f log.*')
            os.system(self.bashDir+'AllrunExceptTheMesh')
            self.mv_dir_and_files(P)
            self.parameters_for_Studies()
    
    def run_PW_first_step(self):
        # Run the first step from nothing 
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        # run the first step
        P  =self.parameters_list[0]
        self.set_parameters_dictionnary(P)
        print(self.para_name+' = '+str(self.parameters[self.para_name]))
        self.print_parameters_for_OF()
        os.system('rm -f log.*')
        os.system(self.bashDir+'AllrunExceptTheMesh')
        #self.cp_dir_and_files(P)
        self.parameters_for_Studies()


################# REV runs functions ##############################
    def run_REV(self):
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        for P in self.parameters_list:
            if self.REV.DoGen:
                self.REV.Cyls = []
            os.system(self.bashDir+'Allclean')
            self.set_parameters_dictionnary(P)
            self.refinement_if_needed(P)
            print(self.para_name+' = '+str(self.parameters[self.para_name]))
            self.eMeshRVE()
            self.Uf = self.U/(1.-self.REV.phi)
            print("mean fluid vel :",self.Uf)
            self.set_parameters_dictionnary(P)
            self.print_parameters_for_OF()
            os.system(self.bashDir+'Allrun1 '+str(self.whichMesh))
            self.printPatch()
            self.forcesFunc('p','U')
            os.system(self.bashDir+'Allrun2 '+str(self.whichMesh))
            self.cp_dir_and_files(P)
            self.Cyls_for_Studies(P)
            self.parameters_for_Studies()
            

    def run_REV_first_step(self):
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir)
        os.system(self.bashDir+'Allclean')
        P  =self.parameters_list[0]
        self.set_parameters_dictionnary(P)
        self.refinement_if_needed(P)
        print(self.para_name+' = '+str(self.parameters[self.para_name]))
        self.eMeshRVE()
        self.Uf = self.U/(1.-self.REV.phi)
        print("mean fluid vel :",self.Uf)
        self.set_parameters_dictionnary(P)
        self.print_parameters_for_OF()
        os.system(self.bashDir+'Allrun1 '+str(self.whichMesh))
        self.printPatch()
        self.forcesFunc('p','U')
        os.system(self.bashDir+'Allrun2 '+str(self.whichMesh))
        self.mv_dir_and_files(P)
        self.parameters_for_Studies()
        self.Cyls_for_Studies(P)

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
        self.nu = parameters['r']*2*parameters['U']/parameters['Re']
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


    
    def update_parameters_dictionnary(self):
        self.parameters = {
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
            "phi"       :self.REV.phi
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
            M = zip(M_x,M_y,M_z)
            name="cercle"+str(numero_du_cercle)
            self.eMesh(M,name,1)
        self.eMeshForSnappy()
    
    def eMeshRVE(self):
        if self.doREV == 1:
            # self.REV.generation()
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
                self.forces.append('\ttype            forces;')
                self.forces.append('\tlibs            (forces);')
                self.forces.append('\twriteControl    timeStep;')
                self.forces.append('\twriteInterval   $EndTime;')
                if self.wrinteAtEnd:
                    self.forces.append('\texecuteControl    timeStep;')
                    self.forces.append('\texecuteInterval   $EndTime;')
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
                self.forces.append('\ttype            forces;')
                self.forces.append('\tlibs            (forces);')
                self.forces.append('\twriteControl    timeStep;')
                self.forces.append('\twriteInterval    $EndTime;')
                if self.wrinteAtEnd:
                    self.forces.append('\texecuteControl    timeStep;')
                    self.forces.append('\texecuteInterval   $EndTime;')
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
        for key,value in self.parameters.iteritems() :
            if type(value) != str:
                parametersOF.append(key+"      "+str(value)+";"+"\n")
        np.savetxt('system/parametersOF',parametersOF, fmt='%s',delimiter=" ")

    def cp_dir_and_files(self,P):
        namedir_simulation = self.namedir+str(P)
        try:
            shutil.rmtree(namedir_simulation)
            os.mkdir(namedir_simulation)
        except:
            os.mkdir(namedir_simulation)
        os.system('cp -r --backup=simple postProcessing/*'+' '+namedir_simulation)
        print os.listdir('postProcessing/.')
        os.system('cp -r --backup=simple log*  '+namedir_simulation)
            # that is for having the name of the last dir
        dirs = os.listdir('.')
        dirs = filter(lambda x : x!='' and x !='.' and x!='0.' , [re.sub("[^0123456789\.]","",dir) for dir in dirs])
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
            shutil.rmtree(namedir_simulation)
            os.mkdir(namedir_simulation)
        except:
            os.mkdir(namedir_simulation)
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
        """print the parameter C++ file for OpenFoam"""
        files = []    
        files.append('parameters = {')
        for key,value in self.parameters.iteritems() :
            if type(value) != str:
                files.append('  "'+key+'":'+str(value)+',')
            else:
                files.append('  "'+key+'":'+'"'+str(value)+'"'+',')
        # files.append('  "para_name":'+'"'+str(self.para)+'"'+',')
        # files.append('  "namedir":'+'"'+str(namedir)+'"')
        files.append('}')
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
        np.savetxt(namedir_simulation+'/Cyls_for_this_study.py',files, fmt='%s',delimiter=" ")
        np.savetxt('Cyls_for_this_study.py',files, fmt='%s',delimiter=" ")
        
    def refinement_if_needed(self,P):
        if self.Auto_refinement[0] == 1 and self.Maillage_de_fond >= self.Auto_refinement[1]:
            self.ER = self.ER + 1
            self.Raffinement_de_surface = self.Raffinement_de_surface + 1
            self.Raffinement_de_region = self.Raffinement_de_region + 1
            self.set_parameters_dictionnary(P)
    
    def cp_first_step_here(self,P):
        '''this function take the constant and last step to the case'''
        namedir_simulation = self.namedir+str(P)+'/'
        # getting the last dir 
        dirs = os.listdir(namedir_simulation)
        dirs = filter(lambda x : x!='' and x !='.' and x!='0.' , [re.sub("[^0123456789\.]","",dir) for dir in dirs])
        dirs =[dire for dire in dirs if dire in os.listdir(namedir_simulation)]
        dirs.sort(key = lambda x : float(x))#range dans l'ordre croissant 
        laststep = dirs[-1]

        laststep_file = namedir_simulation + laststep

        # mv the files here 
        os.system('cp -r '+laststep_file+' '+laststep)
    
    def cp_constant_file_here(self,P):        
        '''this function take the constant dir to the case'''
        namedir_simulation = self.namedir+str(P)+'/'
        constant_file = namedir_simulation + 'constant'
        os.system('rm -rf constant/')

        os.system('cp -rf  '+constant_file+' constant')
