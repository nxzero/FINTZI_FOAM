import re
import os
import math
import sys
import os
import shutil
import re
import numpy as np
from numpy import linalg as LA
# import matplotlib.pyplot as plt
import py_functions.moment_Cox as mc
class Study():
    def __init__(self,namestud,forcestopbot = False,maindir='forces'):
    # def Postprocess(self, namedir):
        self.namestud = namestud
        self.forcestopbot = forcestopbot
        if not namestud == '':
            self.namedir = '../results/'+namestud+'/'
            sys.path.append('./'+self.namedir)
            import parameters_for_this_study as pstud
            from importlib import reload
            reload(pstud)
            self.parameters = pstud.parameters
            sys.path.remove('./'+self.namedir)
            self.maindir = maindir
            self.studies_datas = {}
            self.para_name = self.parameters['para_name']
            self.namesOrig = os.listdir(self.namedir)
            self.names = os.listdir(self.namedir)
            self.names = filter(lambda x : x!='' and x !='.'  , [re.sub("[^0123456789\.]","",dir) for dir in self.names])
            self.names =[dire for dire in self.names if dire in os.listdir(self.namedir)]
            self.names = [dir for dir in set(self.names) if dir in self.namesOrig]
            self.names.sort(key = lambda x : float(x)) #sort par ordre croissant
            self.parameters_list = [float(name) for name in self.names]
            self.name = ''
            self.ori = np.array([6,6,6])
        else:
            self.names = []
            self.initlist()
    

    def get_all_the_datas(self):
        self.get_forces_and_troques()
        self.get_convergences_datas()

    def get_studies_last_step(self):
        for name in self.names :
            self.get_study_last_step(name,self.maindir)
            self.calcul(name)

    def get_study_last_step(self,name,maindir):        
        #definition des parametres
        dirF = self.namedir+str(name)+'/'+maindir+'/'
        dirs = os.listdir(dirF)
        dirs = filter(lambda x : x!='' and x !='.' and x!='0.' , [re.sub("[^0123456789\.]","",dir) for dir in dirs])
        dirs =[dire for dire in dirs if dire in os.listdir(dirF)]
        dirs.sort(key = lambda x : float(x))#range dans l'ordre croissant 
        lastdir = dirs[-1]
        Forces_path = dirF + lastdir +'/force.dat'
        Forces = open(Forces_path,'r').readlines()[-1].split()
        Forces = [float(re.sub("[^0123456789\.e+-]","",force)) for force in Forces]
        self.Forces = {
            'total' : {
                'x' : Forces[1],
                'y' : Forces[2],
                'z' : Forces[3]
            },
            'pressure' :{
                'x' : Forces[4],
                'y' : Forces[5],
                'z' : Forces[6]
            },
            'viscous':{
                'x' : Forces[7],
                'y' : Forces[8],
                'z' : Forces[9]
            }
        }
        Torques_path = dirF + lastdir +'/moment.dat'
        Torques = open(Torques_path,'r').readlines()[-1].split()
        Torques = [float(re.sub("[^0123456789\.e+-]","",Torque)) for Torque in Torques]
        self.Torques = {
            'total' : {
                'x' : Torques[1],
                'y' : Torques[2],
                'z' : Torques[3]
            },
            'pressure' :{
                'x' : Torques[4],
                'y' : Torques[5],
                'z' : Torques[6]
            },
            'viscous':{
                'x' : Torques[7],
                'y' : Torques[8],
                'z' : Torques[9]
            }
        }
        self.studies_datas[name] = {
            'Forces' : self.Forces,
            'Torques' : self.Torques,
            'parameters' : self.parameters,
            'para_value' : float(name),
            'namedir':self.namedir
        }        
    

        

    def calcul(self,name):
        import math
        ############ maj du parametre variable ################
        self.studies_datas[name]['parameters'][self.para_name] = self.studies_datas[name]['para_value']
        ###########Re calculation of the parameters witch depends on the para_value
        self.studies_datas[name]['Re'] =  self.studies_datas[name]['parameters']['Re']
        self.studies_datas[name]['r'] = self.studies_datas[name]['parameters']['r']
        self.studies_datas[name]['nu'] = self.studies_datas[name]['parameters']['r']*2*self.studies_datas[name]['parameters']['U']/self.studies_datas[name]['parameters']['Re']
        self.studies_datas[name]['eta'] = self.studies_datas[name]['nu']*self.studies_datas[name]['parameters']['Rho']
        self.studies_datas[name]['rho'] = self.studies_datas[name]['parameters']['Rho']
        self.studies_datas[name]['length'] = self.studies_datas[name]['parameters']['ksi']*self.studies_datas[name]['parameters']['r']*2
        ################ normalisation de la force :
        self.studies_datas[name]['F_h'] = 6*self.studies_datas[name]['eta']*math.pi*self.studies_datas[name]['parameters']['U']*(3./4.*self.studies_datas[name]['parameters']['r']**2*self.studies_datas[name]['length'])**(1./3.)
        self.studies_datas[name]['F_ad'] = self.studies_datas[name]['Forces']['total']['x']/self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_dp'] = self.studies_datas[name]['Forces']['pressure']['x']/self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_dmu'] = self.studies_datas[name]['Forces']['viscous']['x']/self.studies_datas[name]['F_h']
        ############### Dragcoef #################
        geometrycalparameters = 1./2.*self.studies_datas[name]['parameters']['Rho']*self.studies_datas[name]['parameters']['r']*2*self.studies_datas[name]['parameters']['U']**2*self.studies_datas[name]['parameters']['length']
        self.studies_datas[name]['C_d'] = self.studies_datas[name]['Forces']['total']['x']/geometrycalparameters
        self.studies_datas[name]['C_p'] = self.studies_datas[name]['Forces']['pressure']['x']/geometrycalparameters
        self.studies_datas[name]['C_mu'] = self.studies_datas[name]['Forces']['viscous']['x']/geometrycalparameters
        ###############Force paralle and perp
        thetaRad = self.studies_datas[name]['parameters']['Theta']/360.0 * 2.0*math.pi
        self.studies_datas[name]['F_para'] = self.studies_datas[name]['Forces']['total']['x']*math.cos(thetaRad)+self.studies_datas[name]['Forces']['total']['y']*math.sin(thetaRad)
        self.studies_datas[name]['F_perp'] = - self.studies_datas[name]['Forces']['total']['x']*math.sin(thetaRad)+self.studies_datas[name]['Forces']['total']['y']*math.cos(thetaRad)
        self.studies_datas[name]['F_dpara'] = self.studies_datas[name]['F_para']/self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_dperp'] = self.studies_datas[name]['F_perp']/self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_dd'] = self.studies_datas[name]['Forces']['total']['x'] / self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_dd_PW'] = self.studies_datas[name]['Forces']['total']['z'] / self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_dl'] = self.studies_datas[name]['Forces']['total']['y'] / self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_dl_PW'] = self.studies_datas[name]['Forces']['total']['y'] / self.studies_datas[name]['F_h']
        ############ Adimensionnement du moment par eta U L^2 * 0.25
        self.studies_datas[name]['M_had'] = self.studies_datas[name]['Torques']['total']['z'] / (self.studies_datas[name]['eta']*self.studies_datas[name]['length']**2*self.studies_datas[name]['parameters']['U'])
        self.studies_datas[name]['M_hadP'] = self.studies_datas[name]['Torques']['viscous']['z'] / (self.studies_datas[name]['eta']*self.studies_datas[name]['length']**2*self.studies_datas[name]['parameters']['U'])
        self.studies_datas[name]['M_hadV'] = self.studies_datas[name]['Torques']['pressure']['z'] / (self.studies_datas[name]['eta']*self.studies_datas[name]['length']**2*self.studies_datas[name]['parameters']['U'])
        self.studies_datas[name]['M_had_PW'] = self.studies_datas[name]['Torques']['total']['x'] / (self.studies_datas[name]['eta']*self.studies_datas[name]['length']**2*self.studies_datas[name]['parameters']['U'])
        adim  = self.studies_datas[name]['eta']*self.studies_datas[name]['length']**2*self.studies_datas[name]['parameters']['U']
        self.studies_datas[name]['M_x'] = self.studies_datas[name]['Torques']['total']['x']/ adim
        self.studies_datas[name]['M_y'] = self.studies_datas[name]['Torques']['total']['y']/ adim
        self.studies_datas[name]['M_z'] = self.studies_datas[name]['Torques']['total']['z']/ adim
        #################### Force par unit of lengh ############
        self.studies_datas[name]['fx'] = self.studies_datas[name]['Torques']['total']['x']/self.studies_datas[name]['length']
        # force adimensionne par la eta U L/2
        adimensionneur = self.studies_datas[name]['eta'] * self.studies_datas[name]['parameters']['U'] *self.studies_datas[name]['length'] 
        self.studies_datas[name]['F_d'] = self.studies_datas[name]['Forces']['total']['x'] / adimensionneur
        self.studies_datas[name]['F_d_PW'] = self.studies_datas[name]['Forces']['total']['x'] / adimensionneur
        self.studies_datas[name]['F_l'] = self.studies_datas[name]['Forces']['total']['y'] / adimensionneur   
        self.studies_datas[name]['F_l_PW'] = self.studies_datas[name]['Forces']['total']['y'] / adimensionneur   
        self.studies_datas[name]['F_d2para'] = self.studies_datas[name]['F_para']/adimensionneur
        self.studies_datas[name]['F_d2perp'] = self.studies_datas[name]['F_perp']/adimensionneur
        #### la meme mais pour viscous
        # Calcul de F_d avec thetaU
        self.studies_datas[name]['F_dV'] = self.studies_datas[name]['Forces']['viscous']['x']/ adimensionneur
        self.studies_datas[name]['F_dP'] = self.studies_datas[name]['Forces']['pressure']['x']/ adimensionneur
        self.studies_datas[name]['F_lV'] = self.studies_datas[name]['Forces']['viscous']['y']/ adimensionneur
        self.studies_datas[name]['F_lP'] = self.studies_datas[name]['Forces']['pressure']['y']/ adimensionneur
       
        try: 
            self.studies_datas[name]['parameters']['ThetaU']
            x = self.studies_datas[name]['Forces']['total']['x'] / adimensionneur
            xViscous = self.studies_datas[name]['Forces']['viscous']['x'] / adimensionneur
            yViscous = self.studies_datas[name]['Forces']['viscous']['y'] / adimensionneur
            xPressure = self.studies_datas[name]['Forces']['pressure']['x'] / adimensionneur
            yPressure = self.studies_datas[name]['Forces']['pressure']['y'] / adimensionneur
            x2 = self.studies_datas[name]['Forces']['total']['x'] / self.studies_datas[name]['F_h']
            y = self.studies_datas[name]['Forces']['total']['y'] / adimensionneur
            z = self.studies_datas[name]['Forces']['total']['z'] / adimensionneur
            y2 = self.studies_datas[name]['Forces']['total']['y'] / self.studies_datas[name]['F_h']
            z2 = self.studies_datas[name]['Forces']['total']['z'] / self.studies_datas[name]['F_h']
            ThetaU = self.studies_datas[name]['parameters']['ThetaU']/360.0 * 2.0*math.pi
            Theta = self.studies_datas[name]['parameters']['Theta']
            
            self.studies_datas[name]['F_d']         =  x *         math.cos(ThetaU) + y*           math.sin(ThetaU)
            self.studies_datas[name]['F_d_PW']      =  z *      math.cos(ThetaU) + y*           math.sin(ThetaU)
            self.studies_datas[name]['F_dd_PW']     =  z2 *      math.cos(ThetaU) + y2*           math.sin(ThetaU)
            self.studies_datas[name]['F_dd']        =  x2 *        math.cos(ThetaU) + y2*          math.sin(ThetaU)
            self.studies_datas[name]['F_dV']        =  xViscous *  math.cos(ThetaU) + yViscous*    math.sin(ThetaU)
            self.studies_datas[name]['F_dP']        =  xPressure * math.cos(ThetaU) + yPressure*   math.sin(ThetaU)
            self.studies_datas[name]['F_dl']        =  x2 *       math.sin(ThetaU)  - y2*          math.cos(ThetaU)
            self.studies_datas[name]['F_l_PW']      =  z *     math.sin(ThetaU)  - y*           math.cos(ThetaU)
            self.studies_datas[name]['F_dl_PW']     =  z2 *     math.sin(ThetaU)  - y2*           math.cos(ThetaU)
            self.studies_datas[name]['F_l']         = -  x *        math.sin(ThetaU)  + y*           math.cos(ThetaU)
            self.studies_datas[name]['F_lV']        = - xViscous * math.sin(ThetaU) + yViscous*    math.cos(ThetaU)
            self.studies_datas[name]['F_lP']        = - xPressure* math.sin(ThetaU) + yPressure*   math.cos(ThetaU)
            self.studies_datas[name]['M_had']       =  - self.studies_datas[name]['Torques']['total']['z'] / (self.studies_datas[name]['eta']*self.studies_datas[name]['length']**2*self.studies_datas[name]['parameters']['U'])
        
        except:
            1

        #relative error compared to the theory 


        re = self.studies_datas[name]['Re']
        ksi = self.studies_datas[name]['parameters']['ksi']
        theta = self.studies_datas[name]['parameters']['Theta']
        if theta == 0:
            theta += 1e-3
        if theta == 90:
            theta -= 1e-3
        # self.studies_datas[name]['eF_d'] = np.abs(self.studies_datas[name]['F_d'] - mc.D2(re,ksi,theta))/mc.D2(re,ksi,theta)
        # self.studies_datas[name]['eF_l'] = np.abs(self.studies_datas[name]['F_l'] - mc.L2(re,ksi,theta))/mc.L2(re,ksi,theta)
        # self.studies_datas[name]['eM_had'] = np.abs(self.studies_datas[name]['M_had'] - mc.T_CC(re,ksi,theta))/mc.T_CC(re,ksi,theta)
        self.studies_datas[name]['Re_L'] = re*ksi/2
        if self.studies_datas[name]['F_l']:
            self.studies_datas[name]['eF_d'] = np.abs(self.studies_datas[name]['F_d'] - mc.D2(re,ksi,theta))/self.studies_datas[name]['F_d']
            self.studies_datas[name]['eF_l'] = (np.abs(self.studies_datas[name]['F_l']) - np.abs(mc.L2(re,ksi,theta)))/np.abs(self.studies_datas[name]['F_l'])
            self.studies_datas[name]['eM_had'] = (np.abs(self.studies_datas[name]['M_had']) -np.abs( mc.T_CC(re,ksi,theta)))/np.abs(self.studies_datas[name]['M_had'])
        else:
            self.studies_datas[name]['eF_d'] =0# np.abs(self.studies_datas[name]['F_d'] - mc.D2(re,ksi,theta))/self.studies_datas[name]['F_d']
            self.studies_datas[name]['eF_l'] =0# (np.abs(self.studies_datas[name]['F_l']) - np.abs(mc.L2(re,ksi,theta)))/np.abs(self.studies_datas[name]['F_l'])
            self.studies_datas[name]['eM_had'] =0# (np.abs(self.studies_datas[name]['M_had']) -np.abs( mc.T_CC(re,ksi,theta)))/np.abs(self.studies_datas[name]['M_had'])
            

        # For tri priodic cases
        if self.ori[0] != 6:
            self.calcul_for_tripp(name)

    def calcul_for_tripp(self,name):
        adimensionneur = self.studies_datas[name]['eta'] * self.studies_datas[name]['parameters']['U'] *self.studies_datas[name]['length'] 
        adim  = self.studies_datas[name]['eta']*self.studies_datas[name]['length']**2*self.studies_datas[name]['parameters']['U']
        axis = np.cross(np.array([1,0,0]), self.ori)
        axis /= LA.norm(axis,2)
        thetaRad  = np.arccos(abs(self.ori[0]))
        Torque  = [self.studies_datas[name]['Torques']['total']['x'],self.studies_datas[name]['Torques']['total']['y'],self.studies_datas[name]['Torques']['total']['z']]
        self.studies_datas[name]['M_had'] = np.abs(np.dot(Torque,axis)) / adim

        Fl = math.sqrt(self.studies_datas[name]['Forces']['total']['y']**2+self.studies_datas[name]['Forces']['total']['z']**2)
        self.studies_datas[name]['F_para'] = self.studies_datas[name]['Forces']['total']['x']*math.cos(thetaRad)+Fl*math.sin(thetaRad)
        self.studies_datas[name]['F_perp'] = - self.studies_datas[name]['Forces']['total']['x']*math.sin(thetaRad)+Fl*math.cos(thetaRad)
        self.studies_datas[name]['F_dpara'] = self.studies_datas[name]['F_para']/self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_dperp'] = self.studies_datas[name]['F_perp']/self.studies_datas[name]['F_h']
        self.studies_datas[name]['F_l'] = Fl / adimensionneur

    def get_forces_and_troques(self):
        print(self.namestud,self.para_name,self.names)
        self.get_studies_last_step()
        self.initlist()
        for name in self.names:
            self.makelist(name)
        return self

    def initlist(self):
        self.parameters_list = []
        self.Torques = {
            'total' : {
                'x' : [],
                'y' : [],
                'z' : []
            },
            'pressure' :{
                'x' : [],
                'y' : [],
                'z' : []
            },
            'viscous':{
                'x' : [],
                'y' : [],
                'z' : []
            }
        }
        self.Forces = {
            'total' : {
                'x' : [],
                'y' : [],
                'z' : []
            },
            'pressure' :{
                'x' : [],
                'y' : [],
                'z' : []
            },
            'viscous':{
                'x' : [],
                'y' : [],
                'z' : []
            }
        }
        self.F_ad = []
        self.F_h = []
        self.C  = {
            "d":[],
            "p":[],
            "mu":[]
        }
        self.F_dpara = []
        self.F_d2para = []
        self.F_para = []
        self.F_dperp = []
        self.F_d2perp = []
        self.F_perp = []
        self.M_had = []
        self.M_hadP = []
        self.M_hadV = []
        self.eM_had = []
        self.M_had_PW = []
        self.fx = []
        self.F_d = []
        self.eF_d = []
        self.F_d_PW = []
        self.F_l = []
        self.eF_l = []
        self.F_l_PW = []
        self.F_dd = []
        self.F_dd_PW = []
        self.F_dl = []
        self.F_dl_PW = []
        self.F_dV = []
        self.F_dP = []
        self.F_lV = []
        self.F_lP = []
        self.M_x = []
        self.M_y = []
        self.M_z = []
        self.Re_L = []

    def makelist(self,name):
        # if float(name) == 0:
        self.parameters_list.append(self.studies_datas[name]['para_value'])
        for key in self.studies_datas[name]['Torques']:
            for key2,value2 in self.studies_datas[name]['Torques'][key].items():
                self.Torques[key][key2].append(value2)
        for key in self.studies_datas[name]['Forces']:
            for key2,value2 in self.studies_datas[name]['Forces'][key].items():
                self.Forces[key][key2].append(value2)
        self.F_ad.append(self.studies_datas[name]['F_ad'])
        self.F_h.append(self.studies_datas[name]['F_h'])
        self.C['d'].append(self.studies_datas[name]['C_d'])
        self.C['p'].append(self.studies_datas[name]['C_p'])
        self.C['mu'].append(self.studies_datas[name]['C_mu'])
        self.F_dpara.append(self.studies_datas[name]['F_dpara'])
        self.F_d2para.append(self.studies_datas[name]['F_d2para'])
        self.F_para.append(self.studies_datas[name]['F_para'])
        self.F_dperp.append(self.studies_datas[name]['F_dperp'])
        self.F_d2perp.append(self.studies_datas[name]['F_d2perp'])
        self.F_perp.append(self.studies_datas[name]['F_perp'])
        self.M_had.append(self.studies_datas[name]['M_had'])
        self.M_hadP.append(self.studies_datas[name]['M_hadP'])
        self.M_hadV.append(self.studies_datas[name]['M_hadV'])
        self.eM_had.append(self.studies_datas[name]['eM_had'])
        self.M_had_PW.append(self.studies_datas[name]['M_had_PW'])
        self.fx.append(self.studies_datas[name]['fx'])
        self.F_d.append(self.studies_datas[name]['F_d'])
        self.eF_d.append(self.studies_datas[name]['eF_d'])
        self.F_d_PW.append(self.studies_datas[name]['F_d_PW'])
        self.F_l.append(self.studies_datas[name]['F_l'])
        self.eF_l.append(self.studies_datas[name]['eF_l'])
        self.F_l_PW.append(self.studies_datas[name]['F_l_PW'])
        self.F_dd.append(self.studies_datas[name]['F_dd'])
        self.F_dd_PW.append(self.studies_datas[name]['F_dd_PW'])
        self.F_dl.append(self.studies_datas[name]['F_dl'])
        self.F_dl_PW.append(self.studies_datas[name]['F_dl_PW'])
        ### viscous and pressures
        self.F_dV.append(self.studies_datas[name]['F_dV'])
        self.F_dP.append(self.studies_datas[name]['F_dP'])
        self.F_lV.append(self.studies_datas[name]['F_lV'])
        self.F_lP.append(self.studies_datas[name]['F_lP'])
        self.M_x.append(self.studies_datas[name]['M_x'])
        self.M_y.append(self.studies_datas[name]['M_y'])
        self.M_z.append(self.studies_datas[name]['M_z'])
        self.Re_L.append(self.studies_datas[name]['Re_L'])

    def get_convergences_datas(self):
        import re
        import os
        self.convergences_datas = {}
        for name in self.names :
            self.convergences_datas[name] = {}
            self.convergences_datas[name]['iterations'] = []
            #definition des parametres
            files_path = self.namedir+str(name)+'/logs/'
            files = os.listdir(files_path)
            files.remove('foamLog.awk')
            for file in files:
                self.convergences_datas[name][file] = []
                the_file_path = files_path + file
                datas = open(the_file_path,'r').readlines()
                datas = [data.split() for data in datas]

                for data in datas:
                    self.convergences_datas[name][file].append(data[1])
                    if file == 'clockTime_0':
                        self.convergences_datas[name]['iterations'].append(data[0])

    def get_probes_datas(self):
        import re
        import os
        self.probes_datas = {}
        for name in self.names :
            self.probes_datas[name] = {}
            self.probes_datas[name]['iterations'] = []
            #definition des parametres
            files_path = self.namedir+str(name)+'/probes/0/U'
            files = open(files_path).readlines()
            files = [line.split() for line in files]
            self.n_probe = 0
            while files[self.n_probe][1] != 'Time':
                self.n_probe = self.n_probe + 1
            self.n_probe = self.n_probe - 1
            print(self.n_probe)
            self.probes_datas[name]['U'] = {}
            # for line in files:
            #     self.probes_datas[name]['U']['x'] = float(re.sub("[^0123456789\.e+-]","",line[1]))



            # print(probes_datas[name]['U']['x'])
            # [float(re.sub("[^0123456789\.e+-]","",Torque)) for Torque in Torques]


    def parametric_plot(self,Y,style = '-'):
        import matplotlib.pyplot as plt
        plt.plot()
        plt.plot(self.parameters_list, Y,style)
        plt.xlabel(self.para_name)  # Add an x-label to the axes.

    def parametric_plot_Forces(self,axis,what='total',style = '-'):
        import matplotlib.pyplot as plt
        Y=self.Forces['total'][axis]
        Y1=self.Forces['viscous'][axis]
        Y2=self.Forces['pressure'][axis]
        if what =='total':
            self.parametric_plot(Y,style)
        elif what== 'all':
            self.parametric_plot(Y,style)
            self.parametric_plot(Y1,style)
            self.parametric_plot(Y2,style)
            plt.legend(['total','vicous','pressure'])
        plt.ylabel('Forces hydrodynamique sur X [N]')

    def parametric_plot_Torques(self,axis,what='total',style = '-'):
        import matplotlib.pyplot as plt
        Y=self.Torques['total'][axis]
        Y1=self.Torques['viscous'][axis]
        Y2=self.Torques['pressure'][axis]
        if what =='total':
            self.parametric_plot(Y,style)
        elif what== 'all':
            self.parametric_plot(Y,style)
            self.parametric_plot(Y2,style)
            self.parametric_plot(Y,style)
            plt.legend(['total','vicous','pressure'])
        plt.ylabel("Moment hydrodynamique sur "+axis+" en [Nm]")
        

    def convergences_plot(self,name,file,legend):
        import matplotlib.pyplot as plt
        plt.plot(self.convergences_datas[name]['iterations'], self.convergences_datas[name][file],label=legend)  # Plot some data on the axes.
        plt.xlabel('iteration')  
        plt.legend()

    def convergences_plot_U(self,name):
        import matplotlib.pyplot as plt
        self.convergences_plot(name,'Ux_0','Ux')
        self.convergences_plot(name,'Uy_0','Uy')
        self.convergences_plot(name,'Uz_0','Uz')
        plt.yscale('log')

    def convergences_plot_p(self,name):
        import matplotlib.pyplot as plt
        self.convergences_plot(name,'p_0','p1')
        plt.yscale('log')
        plt.ylabel('Residu de p')

        
    def __add__(self, other):
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] += other.studies_datas[self.name]['Torques'][key][key2]
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] += other.studies_datas[self.name]['Forces'][key][key2]
        return self
        
    def __sub__(self, other):
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] -= other.studies_datas[self.name]['Torques'][key][key2]
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] -= other.studies_datas[self.name]['Forces'][key][key2]
        return self


    def __div__(self, other):
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] /= other
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] /= other
        return self

    def __mul__(self, other):
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] *= other.studies_datas[self.name]['Torques'][key][key2]
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] *= other.studies_datas[self.name]['Forces'][key][key2]
        return self

    def __pow__(self, other):
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] = self.studies_datas[self.name]['Torques'][key][key2]**2
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] = self.studies_datas[self.name]['Forces'][key][key2]**2
        
        for key in self.studies_datas[self.name]:
            if key not in ['Torques','Forces','parameters','para_value','namedir','Re','eta','nu','length']:
                self.studies_datas[self.name][key] = self.studies_datas[self.name][key]**2

        
        
        return self

    def sqrt(self):
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] = math.sqrt(self.studies_datas[self.name]['Torques'][key][key2])
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] = math.sqrt(self.studies_datas[self.name]['Forces'][key][key2])
        for key in self.studies_datas[self.name]:
            if key not in ['Torques','Forces','parameters','para_value','namedir','Re','eta','nu','length']:
                self.studies_datas[self.name][key] = math.sqrt( self.studies_datas[self.name][key])

        return self


    def set_to_0(self,other):        
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] = 0
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] = 0
        for key in other.studies_datas[self.name]:
            if key not in ['Torques','Forces','parameters','para_value','namedir','Re','eta','nu','length']:
                self.studies_datas[self.name][key] = 0
        return self

    def divall(self, other):
        if other == 0 :
            return self
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] /= other
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] /= other

        for key in self.studies_datas[self.name]:
            if key not in ['Torques','Forces','parameters','para_value','namedir','Re','eta','nu','length']:
                self.studies_datas[self.name][key] /= other

        return self

    def addall(self, other):
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] += other.studies_datas[self.name]['Torques'][key][key2]
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] += other.studies_datas[self.name]['Forces'][key][key2]

        for key in other.studies_datas[self.name]:
            if key not in ['Torques','Forces','parameters','para_value','namedir','Re','eta','nu','length']:
                self.studies_datas[self.name][key] += other.studies_datas[self.name][key]

        return self

    def suball(self, other):
        for key in self.studies_datas[self.name]['Torques']:
            for key2 in self.studies_datas[self.name]['Torques'][key]:
                self.studies_datas[self.name]['Torques'][key][key2] -= other.studies_datas[self.name]['Torques'][key][key2]
        for key in self.studies_datas[self.name]['Forces']:
            for key2 in self.studies_datas[self.name]['Forces'][key]:
                self.studies_datas[self.name]['Forces'][key][key2] -= other.studies_datas[self.name]['Forces'][key][key2]

        for key in other.studies_datas[self.name]:
            if key not in ['Torques','Forces','parameters','para_value','namedir','Re','eta','nu','length']:
                self.studies_datas[self.name][key] -= other.studies_datas[self.name][key]

        return self

    