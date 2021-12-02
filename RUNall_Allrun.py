import sys
sys.path.append('../myLibs/py_functions')
from AP_class import *
import sphere as G
from create_Cyls_list import POS,ORI
# entr de l'tude parametrique :


A1 = Analyse_Parametrique('Re',[20],'TP/JLtestSNewBC')
A1.REV.shape                    = 'Sphere'
A1.Re                           = 0.02    
A1.ksi 				            = 2
A1.U                            = 1
# Radius of sheres
A1.r                            = 0.2
A1.Theta 			            = 90
A1.ReMesh                       = 1
A1.MeshForHightRe               = 0

# snappyhexMesh
A1.ER                           = 2
A1.Raffinement_de_surface       = 2
A1.Raffinement_de_region        = 2
A1.GapLevel                     = 3
A1.wakezone                     = 0
A1.nCellsBetweenLevels          = 5
A1.ndc                          = 20 # nombre de cells per diametre
# Main Size of the domain 
A1.HowManyD                     = 2
A1.HowManyL                     = 0
# aspect ratio in the x axis
A1.xblockmeshrelatif            = 1




A1.Relaxcoef                    = 0.85
A1.EndTime 			            = 16000 # End time 
A1.REV.gap                      = A1.r *0.05 #  gap minimum between two spheres  
A1.REV.T 			            =[1,1,1]
A1.add_layers                   = 0 
#tol de simpleFoam
A1.ptol                         = 1e-7
A1.Utol                         = 1e-8
A1.pres                         = 1e-6
A1.Ures                         = 1e-7
# For mergin the patches (do not use)
A1.mPatchs                      = 0
#Dp ini 
A1.DeltaP                       = 1
A1.REV.pres                        = 1000 # res of the eMesh

#########################################
# Size of the domain for generation.py // sinon c'est automatique 
# A1.REV.sizex =1
# A1.REV.sizey =1 
# A1.REV.sizez =1   
# POS = [
#     [9.050713440726041448e-01,5.333047077795801671e-01,7.487155453677005745e-01],
#     [3.162117866548386225e-01,4.972760438690805307e-01,4.155473193240496466e-01],
#     [5.380619805920814347e-01,8.041111619813809952e-01,7.192775968292124400e-01],
#     [9.531627901395731683e-01,1.034170161877492111e+00,3.555096572790191756e-01],
#     [3.405516930198635439e-01,6.600087754068728607e-01,1.039850880336744998e+00],
#     [1.050146345968606676e+00,9.334244793438388754e-01,9.572578165705754039e-01],
#     [4.893136914761649914e-01,1.050823348461230600e+00,3.078914081055839080e-01],
#     [5.137868092184737501e-01,3.317077099593062628e-01,8.080484208243485789e-01],
#     [8.439807897294464567e-01,5.304934380124679549e-01,3.405162143138772968e-01],
#     [1.081704566454322558e+00,3.398203403006015422e-01,1.065666790397836960e+00]
# ]

# for pos,ip in zip(POS,range(len(POS))):
#     pos = [pos[0]-0.2,pos[1]-0.2,pos[2]-0.2]
#     A1.REV.Spheres.append([G.SphereObj(pos,0.2),str(ip)])
#     A1.REV.RealSpheres.append([G.SphereObj(pos,0.2),str(ip)])
#########################################
# generation auto

A1.REV.noc                      = len(A1.REV.Spheres)
A1.REV.noc                      = 5
A1.REV.DoGen                    = 1
A1.REV.addPerio                 = 1 #add sheres perio 

A1.REV.DoLoop                   = 1 # shift des spheres tangents 
A1.REV.Rmin                     = 0.9
A1.REV.Rmax                     = 1.1

A1.REV.DoTheRest                = 0 # To do for sphere

A1.REV.DeleteThem               = 0 # To do for sphere 

A1.run_REV()

