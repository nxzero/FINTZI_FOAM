from py_functions.AP_class import *
import math
# parameteric analisys on the Reynolds number
ksi = 2
phi = 0.30
noc = 0
A1 = Analyse_Parametrique('Re',[10,20],f'Deletme/ksi_{ksi}/Phi_{phi}/') 
A1.ksi                          = ksi # mean aspect ratio
A1.HowManyD                     = 4 # size of the domain relative to the
# diameter of the cyl
A1.ndc                          = 10 # number of cells per diameters default 30
A1.ER                           = 3 # refinement of the edges
A1.Relaxcoef                    = 0.95 #Relaxation coefficient
A1.ReMesh                       = 1  
A1.MeshForHightRe               = 0
A1.U                            = 1
A1.Raffinement_de_surface       = 3
A1.Raffinement_de_region        = 3
A1.wakezone                     = 0
A1.nCellsBetweenLevels          = 5 #20
A1.HowManyL                     = 0
A1.xblockmeshrelatif            = 2
A1.EndTime                      = 10
A1.REV.gap                      = A1.r *0.005
A1.add_layers                   = 0
A1.wrinteAtEnd                  = 1
A1.ptol                         = 1e-7
A1.Utol                         = 1e-8
A1.pres                         = 1e-6
A1.Ures                         = 1e-7
A1.UtolBC                       = 0.001
A1.mPatchs                      = 0
A1.REV.noc                      = 1 #number of cyl minus 1
A1.REV.T                        =[1,1,1] # orientation tensor
A1.REV.phi                      = phi
A1.REV.DoREV                    = 1
A1.REV.DoGen                    = 1   # Do generation if True.
A1.REV.DoLoop                   = 1   # Do the loop where the 
#cylinders are shifted. 
A1.REV.DoTheRest                = 1   # Rotate independently 
#the cylinders tengent to a boundary.
A1.REV.DeleteThem               = 0   # Delet remaining cylinder
#if True. 
# Criterion of max and min angles of tangency. 
A1.REV.min2A                    = 20. /180. * math.pi  
A1.REV.minA                     = 10. /180. * math.pi  
A1.REV.min2B                    = 60. /180. * math.pi  
A1.REV.minB                     = 80. /180. * math.pi    

A1.run_a_lots_mesh()
A1.run_parametric() # Run parametrical sim from this mesh 

