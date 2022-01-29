import math
import numpy as np
from numpy import linalg as LA
from py_functions.cylinder import *
from py_functions.sphere import *

class Six_pack():
    def __init__(self):
        self.sizex = 1
        self.sizey =     self.sizez = 1
        self.Matiere = np.array([self.sizex,self.sizey*0.1,self.sizez*0.1])
        self.pRef = np.array([self.sizex/2.,self.sizey/2.,self.sizez - 1e-5])
        self.noc = 10
        self.Cyls = []
        self.Spheres = []
        self.gap = 0
        self.RealCyls = []
        self.RealSpheres = []
        self.rmean = 0.1
        self.lengthmean = 1.
        self.X,self.Y,self.Z = list(),list(),list()

        ########## Physicals parameters ##################
        self.phi = 0.2
        self.T11 = 1./3.
        self.T22 = 1./3.
        self.T33 = 1./3.
        self.T = [self.T11,self.T22,self.T33]
        self.minA = 9./360. * 2.*math.pi
        self.minB = 81./360. * 2.*math.pi
        self.min2A = 15. /360. * 2.*math.pi
        self.min2B = 72. /360. * 2.*math.pi
        self.gapbetweenplanes = 0.005
        self.shape = 'Cyls'
        self.DoLoop = 1
        self.DoGen = 1
        self.DoTheRest =1 
        self.cylsBlocked = []
        self.DeleteThem = 0
        self.CylsMod = set()
        self.CylNoMod = set()
        self.addPerio   =1

    def run(self):        
        self.print_parameters()
        if self.DoGen:
            self.generation()
        else:
            if self.addPerio:
                self.AddselfPerios()
        self.generation_des_planKey() 
        self.Donne_post_generations()
        if self.DoLoop:
            self.DoLoopMethode()
        if self.DoTheRest:
            self.DoTheRestMethode(self.min2A,self.min2B)
        if self.DeleteThem:
            self.DeletMethode()
        self.Donne_post_generations()
        self.generation_des_planKey()
        self.generation_des_lists()
        self.correction_des_lists()   
        self.generation_des_points()

    def DoLoopMethode(self):
        print('##############################################################################')
        print('Shift For : ',self.minB/math.pi*180., 'and ',self.minA/math.pi*180.)
        i = 0
        while self.IsItTangent(self.minA,self.minB):# or self.IsItTooCloseFromTheEdges() == 1:
            i = i + 1
            print(str(i)+" try...")
            if LA.norm(self.Shift,2) >= self.Cylten.L and LA.norm(self.Shift,2):
                self.Shift = self.Shift/LA.norm(self.Shift,2) *self.Cylten.L*0.5
            if LA.norm(self.Shift,2) <= 1e-5:#sometime it bug
                self.Shift = self.Shift/LA.norm(self.Shift,2) *self.Cylten.L*0.5
            self.ShiftAll()
            self.generation_des_planKey()


    def DoTheRestMethode(self,min2A,min2B):
        print('##############################################################################')
        print('Changes To : ',min2B/math.pi*180., 'and ',min2A/math.pi*180.)
        while self.IsItTangent(min2A,min2B):
            if LA.norm(self.Shift,2) <= 1e-5 and LA.norm(self.Shift,2) != 0:#sometime it bug
                self.Cylten.shift = self.Shift/LA.norm(self.Shift,2) *self.Cylten.L*0.25
            if abs(self.Cylten.theta) > min2B:
                self.ChangeCyltenPos(min2A,min2B)
                if self.TryAgain:
                    self.Cylten.shift *= -1
                    self.ChangeCyltenPos(min2A,min2B)
                    self.Cylten.shift *= -1
                    if self.TryAgain:
                        self.ChangeCyltenOri(min2A,min2B)
                        if self.TryAgain:
                            self.ChangeCyltenOri(-min2A,-min2B)
                            if self.TryAgain:
                                thetas1 = np.linspace(-min2B,-min2A,20)
                                thetas2 = np.linspace(min2B,min2A,20)
                                thetas =thetas1.tolist() + thetas2.tolist()
                                for theta in thetas:
                                    if self.TryAgain:
                                        self.ChangeCyltenOri(min2A,theta)
            if abs(self.Cylten.theta) < min2A:

                self.ChangeCyltenOri(min2A,min2B)
                if self.TryAgain:
                    self.ChangeCyltenOri(-min2A,-min2B)
                    if self.TryAgain:
                        self.ChangeCyltenPos(min2A,min2B)
                        if self.TryAgain:
                            self.Cylten.shift *= -1
                            self.ChangeCyltenPos(min2A,min2B)
                            self.Cylten.shift *= -1
                            if self.TryAgain:
                                thetas2 = np.linspace(min2A,min2B,20)
                                thetas1 = np.linspace(-min2A,-min2B,20)
                                thetas =thetas2.tolist()+thetas1.tolist()  
                                for theta in thetas:
                                    if self.TryAgain:
                                        self.ChangeCyltenOri(theta,min2B)

            if self.TryAgain: 
                self.cylsBlocked.extend([cyll[1] for cyll in self.cylsss])
                print ([cyl[1] for cyl in self.cylsss],'CANNOT BE MOVED !!!')
                self.CylNoMod.add(self.cylsss[0][1])


            self.generation_des_planKey()
        Mod  =  len(self.CylsMod)
        NoMod  =  len(self.CylNoMod)
        print(str(Mod)+" Cylinders have been modified and "+str(NoMod)+" cannot...")

    def DeletMethode(self):
        print('##############################################################################')
        print('Deleting...:')
        self.cylsBlocked = []
        i = 0
        while self.IsItTangent(self.min2A,self.min2B):
            self.DeletCyl()
            self.generation_des_planKey()
            i+=1
        print(str(i)+" Cylinders have been deleted...")
        self.noc -= i
    
    def print_parameters(self):
        if self.shape == 'Cyls':
            print( "radius=",self.rmean)
            print( "length=",self.lengthmean)
            print( "size=",self.sizex,self.sizey,self.sizez)
            print( "number of cylinders=",self.noc )
            print( "gap between them=",self.gap )
            print( "Orientation tensor=",self.T)
            print( "gap between planes and cyls=",self.gapbetweenplanes)
        elif self.shape == 'Sphere':
            print( "radius=",self.rmean)
            print( "size=",self.sizex,self.sizey,self.sizez)
            print( "number of sphere=",self.noc )
            print( "gap between them=",self.gap )
            print( "gap between planes and cyls=",self.gapbetweenplanes)

    def generation(self):
        # self.Cyls = []
        if self.shape == 'Cyls':
            self.generationC()
        if self.shape == 'Sphere':
            self.generationS()

    def generationC(self):
        idcyl = 0
        while len(self.RealCyls)<=self.noc:
            phi = np.arccos(1.-2.*np.random.rand())
            the = np.random.rand() *math.pi*2.
            u = np.array([
                np.sin(phi)*np.cos(the)*self.T[0],
                np.sin(phi)*np.sin(the)*self.T[1],
                np.cos(phi)*self.T[2]
            ])
            u = u/LA.norm(u,2)

            p = np.array([np.random.rand()*self.sizex,np.random.rand()*self.sizey,np.random.rand()*self.sizez])
            Cyl =  Cylinder(p,self.rmean,self.lengthmean,u)
            Cyl.gap = self.gap
            Cyl.colisionCheck(self.Cyls)
            if Cyl.notok == 0:

                # create periodic neibough
                self.createCylperio(Cyl)


                Cylstemp = self.Cyls[:]
                Cylstemp.append([Cyl,idcyl])
                Cylstemp.extend(self.CylsPerios)
                for CylPerioname in self.CylsPerios:
                    CylPerio = CylPerioname[0]
                    CylPerio.colisionCheck(Cylstemp)
                    if CylPerio.notok == 1:
                        Cyl.notok = 1
                        print("collide its self")
                        break
                if Cyl.notok == 0:
                    self.Cyls.append([Cyl,str(idcyl)])
                    self.RealCyls.append([Cyl,str(idcyl)])
                    for CylPerioname in self.CylsPerios:
                        CylPerio = CylPerioname[0]
                        name = CylPerioname[1]
                        if CylPerio.OP[0] < self.sizex+CylPerio.L+CylPerio.marge and CylPerio.OP[0] > -CylPerio.L-CylPerio.marge:
                            if CylPerio.OP[1] < self.sizey+CylPerio.L+CylPerio.marge and CylPerio.OP[1] > -CylPerio.L-CylPerio.marge:
                                if CylPerio.OP[2] < self.sizez+CylPerio.L+CylPerio.marge and CylPerio.OP[2] > -CylPerio.L-CylPerio.marge:
                                    idcylperio = str(idcyl) +'_'+ name
                                    self.Cyls.append([CylPerio,idcylperio])
                    idcyl = idcyl+1
                    print(Cyl, Cyl.e1,len(self.Cyls),idcyl)

    def generationS(self):
        idSphere = 0
        while len(self.RealSpheres)<=self.noc:
            p = np.array([np.random.rand()*self.sizex,np.random.rand()*self.sizey,np.random.rand()*self.sizez])
            Sphere =  SphereObj(p,self.rmean)
            Sphere.gap = self.gap
            Sphere.colisionCheck(self.Spheres)
            if Sphere.notok == 0:

                # create periodic neibough
                self.createSphereperio(Sphere)


                Spherestemp = self.Spheres[:]
                Spherestemp.append([Sphere,idSphere])
                Spherestemp.extend(self.SpheresPerios)
                for SpherePerioname in self.SpheresPerios:
                    SpherePerio = SpherePerioname[0]
                    SpherePerio.colisionCheck(Spherestemp)
                    if SpherePerio.notok == 1:
                        Sphere.notok = 1
                        print(Sphere.notok)
                        break
                if Sphere.notok == 0:
                    self.Spheres.append([Sphere,str(idSphere)])
                    self.RealSpheres.append([Sphere,str(idSphere)])
                    for SpherePerioname in self.SpheresPerios:
                        SpherePerio = SpherePerioname[0]
                        name = SpherePerioname[1]
                        if SpherePerio.OP[0] < self.sizex+SpherePerio.r+SpherePerio.marge and SpherePerio.OP[0] > -SpherePerio.r-SpherePerio.marge:
                            if SpherePerio.OP[1] < self.sizey+SpherePerio.r+SpherePerio.marge and SpherePerio.OP[1] > -SpherePerio.r-SpherePerio.marge:
                                if SpherePerio.OP[2] < self.sizez+SpherePerio.r+SpherePerio.marge and SpherePerio.OP[2] > -SpherePerio.r-SpherePerio.marge:
                                    idSphereperio = str(idSphere) +'_'+ name
                                    self.Spheres.append([SpherePerio,idSphereperio])
                    idSphere = idSphere+1
                    print(Sphere, Sphere.OP,len(self.Spheres),idSphere)
    #fonction a utiliser manuellement si besoin
    def AddselfPerios(self):
        print(self.sizex,self.sizey,self.sizez)
        for Cyl in self.Cyls[:]:
            idcyl = Cyl[1]
            Cyl = Cyl[0]
            self.createCylperio(Cyl)
            for CylPerioname in self.CylsPerios:
                CylPerio = CylPerioname[0]
                name = CylPerioname[1]
                # print CylPerio.OP.all() == 0
                if CylPerio.OP[0] < self.sizex+CylPerio.L+CylPerio.marge and CylPerio.OP[0] > -CylPerio.L-CylPerio.marge:
                    if CylPerio.OP[1] < self.sizey+CylPerio.L+CylPerio.marge and CylPerio.OP[1] > -CylPerio.L-CylPerio.marge:
                        if CylPerio.OP[2] < self.sizez+CylPerio.L+CylPerio.marge and CylPerio.OP[2] > -CylPerio.L-CylPerio.marge:
                            idcylperio = str(idcyl) +'_'+ name
                            self.Cyls.append([CylPerio,idcylperio])
                            print(CylPerio, CylPerio.e1,len(self.Cyls),idcylperio)
        for Sphere in self.Spheres[:]:
            idSphere = Sphere[1]
            Sphere = Sphere[0]
            self.createSphereperio(Sphere)
            for SpherePerioname in self.SpheresPerios:
                SpherePerio = SpherePerioname[0]
                name = SpherePerioname[1]
                # print SpherePerio.OP.all() == 0
                if SpherePerio.OP[0] < self.sizex+SpherePerio.r+SpherePerio.marge and SpherePerio.OP[0] > -SpherePerio.r-SpherePerio.marge:
                    if SpherePerio.OP[1] < self.sizey+SpherePerio.r+SpherePerio.marge and SpherePerio.OP[1] > -SpherePerio.r-SpherePerio.marge:
                        if SpherePerio.OP[2] < self.sizez+SpherePerio.r+SpherePerio.marge and SpherePerio.OP[2] > -SpherePerio.r-SpherePerio.marge:
                            idSphereperio = str(idSphere) +'_'+ name
                            self.Spheres.append([SpherePerio,idSphereperio])
                            print(SpherePerio, SpherePerio.OP,len(self.Spheres),idSphereperio)

    def generation_des_planKey(self):
        if self.shape =='Cyls':
            for Cyl in self.Cyls:
                Cyl[0].calcul_planKey(self.sizex,self.sizey,self.sizez)
        elif self.shape =='Sphere':
            for Sphere in self.Spheres:
                Sphere[0].calcul_planKey(self.sizex,self.sizey,self.sizez)
        else:
            print('!!!!! Wrong Analyse_parameteric.REV.shape argument choose between "Cyls" or "Sphere" !!!!!')
            
    def generation_des_lists(self):
        if self.shape =='Cyls':
            for Cyl in self.Cyls:
                Cyl[0].calcul_list()
        elif self.shape =='Sphere':
            for Sphere in self.Spheres:
                Sphere[0].calcul_list()

    def correction_des_lists(self):
        if self.shape =='Cyls':
            for Cyl in self.Cyls:
                Cyl[0].rearrangeList()
        if self.shape =='Sphere':
            for Sphere in self.Spheres:
                Sphere[0].rearrangeList()

    def generation_des_points(self):
        self.Matiere = np.array([self.sizex*0.1,self.sizey*0.1,self.sizez*0.1])
        self.pRef = np.array([self.sizex/2.,self.sizey/2.,self.sizez - 1e-5])
        self.dans_un_cylindre = 1
        self.dans_un_cylindre2 = 1
        self.dans_une_sphere = 1
        self.dans_une_sphere2 = 1
        if self.shape == 'Cyls':
            print("generation des points...")
            # re define matiere if needed
            while self.dans_un_cylindre == 1 or self.dans_un_cylindre2 == 1:
                self.dans_un_cylindre = 0
                self.dans_un_cylindre2 = 0
                for Cyl in [Cyls[0] for Cyls in self.Cyls]:
                    self.dist = Cyl.OP - self.Matiere
                    self.dist2 = Cyl.OP - self.pRef
                    self.distR = np.array([np.dot(self.dist,Cyl.e2),np.dot(self.dist,Cyl.e3)])
                    self.distR2 = np.array([np.dot(self.dist2,Cyl.e2),np.dot(self.dist2,Cyl.e3)])
                    self.distR = LA.norm(self.distR,2)
                    self.distR2 = LA.norm(self.distR2,2)
                    if self.distR < (Cyl.r*1.1) :#and LA.norm(self.dist ,2)< (Cyl.L*1.1):
                        self.dans_un_cylindre = 1  
                    elif self.distR2 < (Cyl.r*1.1):#and LA.norm(self.dist2,2) < (Cyl.L+self.gap):
                        self.dans_un_cylindre2 = 1

                if self.dans_un_cylindre == 1:
                    self.Matiere = np.array([np.random.rand()*self.sizex,np.random.rand()*self.sizey,np.random.rand()*self.sizez])

                if self.dans_un_cylindre2 == 1:
                    randp = np.array([np.random.rand()*self.sizex,np.random.rand()*self.sizey,np.random.rand()*self.sizez])
                    self.pRef = np.array([self.sizex - 1e-5,randp[1],randp[2]])
        if self.shape == 'Sphere':  
            print("generation des points...")
            # re define matiere if needed
            while self.dans_une_sphere== 1 or self.dans_une_sphere2 == 1:
                self.dans_une_sphere = 0
                self.dans_une_sphere2 = 0
                for Sphere in [Spheres[0] for Spheres in self.Spheres]:
                    self.dist = Sphere.OP - self.Matiere
                    self.dist2 = Sphere.OP - self.pRef
                    self.distR = LA.norm(self.dist,2)
                    self.distR2 = LA.norm(self.dist2,2)
                    if self.distR < (Sphere.r*1.1) :#and LA.norm(self.dist ,2)< (Sphere.L+self.gap):
                        self.dans_une_sphere = 1
                    elif self.distR2 < (Sphere.r*1.1):#and LA.norm(self.dist2,2) < (Sphere.L+self.gap):
                        self.dans_une_sphere2 = 1

                if self.dans_une_sphere == 1:
                    self.Matiere = np.array([np.random.rand()*self.sizex,np.random.rand()*self.sizey,np.random.rand()*self.sizez])

                if self.dans_une_sphere2 == 1:
                    randp = np.array([np.random.rand()*self.sizex,np.random.rand()*self.sizey,np.random.rand()*self.sizez])
                    self.pRef = np.array([self.sizex - 1e-5,randp[1],randp[2]])

            print('Matiere= ', self.Matiere)
            print('PRef = ', self.pRef)

    def IsItTangent(self,minA,minB):
        """The goal is to shift all cylinders if one is tangent to the sides"""
        CylsTemp = self.Cyls[:]
        for cylip in self.cylsBlocked:
            index = [cyl2[1] for cyl2 in CylsTemp].index(cylip)
            del CylsTemp[index]
        if self.shape == 'Cyls':
            for Cyl in CylsTemp:
                idCyl = Cyl[1]
                Cyl = Cyl[0]
                Cyl.PC = Cyl.e1 * Cyl.length/2
                Cyl.OC1 = Cyl.OP + Cyl.PC
                Cyl.OC2 = Cyl.OP - Cyl.PC
                for key in Cyl.plansKeys :
                    if key not in ['Ixt','Iyt','Izt']:
                        continue
                    self.ItIsTooTangent = 0 
                    self.ItIsTangentButItIsFine = 0
                    self.TheEdgesIsnotOk = 0
                    if key in Cyl.plansKeys:
                        n = Cyl.plans[key][0]
                        Ip = Cyl.plans[key][1]
                        # distance to the plan :
                        Dist1 = Ip - Cyl.OC1  
                        Dist2 = Ip - Cyl.OC2  
                        # distance normal
                        Dist1 = np.dot(Dist1,n)
                        Dist2 = np.dot(Dist2,n)
                        # Angle d inclinaison par rapport au plan 
                        beta = np.arccos(np.dot(n,Cyl.e1))
                        theta = math.pi/2. - beta
                        #tengent or not really
                        if abs(theta) <= minA:
                            self.ItIsTooTangent = 1 
                        if abs(Dist1) <= Cyl.r*0.96 and abs(Dist2) <= Cyl.r*0.96:
                            self.ItIsTangentButItIsFine = 1
                        if abs(Dist1) >= (Cyl.r +self.gapbetweenplanes) and abs(Dist2) >= (Cyl.r +self.gapbetweenplanes):
                            self.ItIsTangentButItIsFine = 1

                        # edges stuff
                        # if abs(Dist1) >= (Cyl.r)*np.cos(thetaG)*0.95 and abs(Dist1) <= (Cyl.r)*np.cos(theta):
                        #     self.TheEdgesIsnotOk = 1
                        # if abs(Dist2) >= (Cyl.r)*np.cos(theta)*0.95 and abs(Dist2) <= (Cyl.r)*np.cos(theta):
                        #     self.TheEdgesIsnotOk = 1
                        Cyl.theta = theta
                        Cyl.axis = np.cross(Cyl.e1,n)/LA.norm(np.cross(Cyl.e1,n))
                        Cyl.Dist1 = Dist1
                        Cyl.Dist2 = Dist2
                        Cyl.id = idCyl
                        Cyl.shift = (Cyl.r*2*np.cos(theta)+self.gapbetweenplanes) * n
                        #The plate is tangent 
                        if ('H'+key in Cyl.plansKeys) or ('B'+key in Cyl.plansKeys):
                            if abs(theta) >= minB:
                                print('The Cylinder ,',idCyl,'planz is tengent to the plan',key,'with an angle = ',theta*360./(2.*math.pi))
                                self.Shift = Cyl.shift
      
                                self.Cylten = Cyl
                                return 1
                        # shifting
                        # if self.TheEdgesIsnotOk == 1 and Dist1*Dist2 > 0:
                        #     distmin = min(Dist1,Dist2) 
                        #     print('The Cylinder ,',idCyl,'edges is tengent t the plan',key,'with Dist = ',distmin)
                        #     self.Cylten = Cyl
                        #     self.Shift =  (Cyl.r*np.cos(theta)+self.gap) * n
                        #     return 1

                        if self.ItIsTooTangent == 1 and self.ItIsTangentButItIsFine == 0:
                            distmin = min(Dist1,Dist2)
                            distmax = max(Dist1,Dist2)
                            print( 'The Cylinder ',idCyl,'is tangent to the',key,'plan, Theta = ',theta*360./(2.*math.pi))
                            if theta <= Cyl.AngleChill and (Dist1 + Dist2)/2 >= 0:
                                self.Shift =  (Dist1 + Dist2)/2 * n
                                Cyl.shift = self.Shift
                            elif distmax >= 0 and distmin <= 0:
                                self.Shift =  (distmax + Cyl.r) * n
                                Cyl.shift = self.Shift
                            elif distmax <= 0 and distmin <= 0:
                                self.Shift =  (Cyl.r + distmin) * n
                                Cyl.shift = self.Shift
                            elif distmax >= 0 and distmin >= 0:
                                self.Shift =  (Cyl.r + distmax) * n
                                Cyl.shift = self.Shift
                            self.Cylten = Cyl

                            # I want to shift all by some values
                            return 1
                
            return 0
        if self.shape == 'Sphere':
            """The goal is to shift all cylinders if one is tangent to the sides"""
            for Sphere in self.Spheres:
                idSphere = Sphere[1]
                Sphere = Sphere[0]
                for key in Sphere.plansKeys :
                    if key not in ['Ixt','Iyt','Izt']:
                        continue
                    self.ItIsTooTangent = 0 
                    if Sphere.pts[key]['X'] != [] :
                        n = Sphere.plans[key][0]
                        Ip = Sphere.plans[key][1]
                        # distance to the plan :
                        Dist1 = Ip - Sphere.OP  
                        # distance normal
                        Dist1 = np.dot(Dist1,n)
                        # Angle d inclinaison par rapport au plan 
                        if abs(Dist1) >= Sphere.r*0.85 and abs(Dist1) <= Sphere.r*1.01:
                            self.ItIsTooTangent = 1
                        if self.ItIsTooTangent == 1:
                            print( 'The Sphere ',idSphere,'is tangent to the',key,'plan')
                            self.Shift =  (0.1 * Sphere.r) * n

                            self.Cylten = Sphere

                            # I want to shift all by some values
                            return 1
            return 0

    def ChangeCyltenOri(self,min2A,min2B):
        self.TryAgain = 0
        ID= self.Cylten.id.split('_')[0]
        ids = [ids[1] for ids in self.Cyls if ids[1].split('_')[0] == ID]
        self.cylsss = []
        NOOK = 0 
        for cyl in self.Cyls:
            
            if cyl[1] in ids:
                cyl[0].theta =self.Cylten.theta 
                cyl[0].axis  =self.Cylten.axis 
                cyl[0].Dist1 =self.Cylten.Dist1 
                cyl[0].Dist2 =self.Cylten.Dist2 
                cyl[0].shift =self.Cylten.shift 
                if abs(cyl[0].theta) > abs(min2B):
                    thetamin = min2B - 0.0001
                    if min2B < 0:
                        thetamin = -min2B + 0.0001
                else:
                    thetamin = min2A + 0.0001
                    if min2A < 0:
                        thetamin = - min2A - 0.0001
                rot = thetamin - cyl[0].theta
                axisOfrot = cyl[0].axis
                e1  = cyl[0].e1 * np.cos(rot) + (np.cross(axisOfrot,cyl[0].e1))*np.sin(rot)+axisOfrot*(np.  dot(axisOfrot,cyl[0].e1))*(1-np.cos(rot))
                cyltemp = Cylinder(cyl[0].OP,cyl[0].r, cyl[0].length, e1)
                CylsCopy = self.Cyls[:]
                index = [cyl2[1] for cyl2 in self.Cyls].index(cyl[1])
                del CylsCopy[index]
                cyltemp.gap = self.gap
                cyltemp.colisionCheck(CylsCopy)
                if cyltemp.notok:
                    NOOK = 1
                    cylnotok = cyl[1]
                self.cylsss.append([cyltemp,cyl[1]])
                
        if  NOOK == 0 and not self.cylsss[0][1].split('_').count('Mo') > 2:
            for cyl in self.cylsss:
                index = [cyl2[1] for cyl2 in self.Cyls].index(cyl[1])
                del self.Cyls[index]
                cyl[1] += '_Mo'
                cyl[0].lastShift = self.Cylten.shift
                self.Cyls.append(cyl)
            self.CylsMod.add(self.cylsss[0][1])    
            print ([cyl[1] for cyl in self.cylsss],'Has been rotated to ',thetamin/math.pi *180. )
        else:
            self.TryAgain = 1

    


    def ChangeCyltenPos(self,min2A,min2B):
        self.TryAgain = 0
        min2A += 0.001
        min2B -= 0.001
        ID= self.Cylten.id.split('_')[0]
        ids = [ids[1] for ids in self.Cyls if ids[1].split('_')[0] == ID]
        self.cylsss = []
        NOOK = 0 
        for cyl in self.Cyls:
            if cyl[1] in ids:
                cyl[0].theta =self.Cylten.theta 
                cyl[0].axis  =self.Cylten.axis 
                cyl[0].Dist1 =self.Cylten.Dist1 
                cyl[0].Dist2 =self.Cylten.Dist2 
                cyl[0].shift =self.Cylten.shift 
                
                # if cyl[1].split('_').count('Mo') >= 1 and cyl[0].shift = - cyl[0].lastShift:
                    # cyl[0].shift = + cyl[0].lastShift
                OP = cyl[0].OP + cyl[0].shift
                cyltemp = Cylinder(OP,cyl[0].r, cyl[0].length, cyl[0].e1)
                cyltemp.gap = self.gap
                CylsCopy = self.Cyls[:]
                index = [cyl2[1] for cyl2 in self.Cyls].index(cyl[1])
                del CylsCopy[index]
                cyltemp.colisionCheck(CylsCopy)
                if cyltemp.notok:
                    NOOK = 1
                    cylnotok = cyl[1]
                self.cylsss.append([cyltemp,cyl[1]])


        if  NOOK == 0 and not self.cylsss[0][1].split('_').count('Mo') > 2:
            for cyl in self.cylsss:
                index = [cyl2[1] for cyl2 in self.Cyls].index(cyl[1])
                del self.Cyls[index]
                cyl[1] += '_Mo'
                cyl[0].lastShift = self.Cylten.shift
                self.Cyls.append(cyl)    
            self.CylsMod.add(self.cylsss[0][1]) 
            print ([cyl[1] for cyl in self.cylsss],'Has been translate DL = ',self.Cylten.shift )
        else:
            self.TryAgain = 1
            







    def DeletCyl(self):
        ID= self.Cylten.id.split('_')[0]
        ids = [ids[1] for ids in self.Cyls if ids[1].split('_')[0] == ID]
        self.Cyls = [cyl for cyl in self.Cyls if not cyl[1] in ids]
        self.RealCyls = [cyl for cyl in self.RealCyls if not cyl[1] in [ID]]


    def IsItTooCloseFromTheEdges(self):
        ########## define all the 12 verticies ################
        orgib = np.array([0 , 0 , 0 ])
        EdgesVec1 = np.array([1, 0, 0])
        EdgesVec2 = np.array([0, 1, 0])
        EdgesVec3 = np.array([0, 0, 1])
        ########## orig top ##################
        origt = np.array([self.sizex,self.sizey,self.sizez])
        origtz = np.array([self.sizex,self.sizey,0])
        origty = np.array([self.sizex,0,self.sizez])
        origtx = np.array([0,self.sizey,self.sizez])
        points = [orgib,origt,origtx,origty,origtz]
        Directions = [EdgesVec1,EdgesVec2,EdgesVec3]

        if self.shape == 'Cyls':
            for Cyl in self.Cyls:
                idcyl = Cyl[1]
                Cyl = Cyl[0]
                for point in points:
                    for direction in Directions:
                        dist = Cyl.OP - point 
                        ### dist between point and center ###
                    
                        nz = np.cross(dist,direction)
                        nplan = np.cross(direction,nz)
                        nplan = nplan / LA.norm(nplan)
                        dist3 = abs(np.dot(nplan,dist)) 
                        ### distance between axis ###
                        nA = np.cross(Cyl.e1,direction)
                        distA = np.dot(dist,nA)
    
                        if abs(distA) > 0.8*Cyl.r and abs(distA) < 1.1*Cyl.r and dist3 < Cyl.L:
                            print('Cylinder :',idcyl,'too colse from an edges',abs(distA),'dist3',dist3)
                            self.Shift =  0.1*Cyl.r * nA
                            print( nA)
                            self.Cylten = Cyl
                            return 1
            return 0 
        if self.shape == 'Sphere':
            for Sphere in self.Spheres:
                idSphere = Sphere[1]
                Sphere = Sphere[0]
                for point in points:
                    for direction in Directions:
                        dist = Sphere.OP - point 
                        ### dist between point and center ###

                        nz = np.cross(dist,direction)
                        nplan = np.cross(direction,nz)
                        try:
                            nplan = nplan / LA.norm(nplan)
                        except:
                            nplan = np.array([0,0,0])
                        dist3 = abs(np.dot(nplan,dist)) 
                        ### distance between axis ###

                        if abs(dist3) > 0.9*Sphere.r and abs(dist3) < 1.1*Sphere.r:
                            print('Sphere :',idSphere,'too colse from an edges dist = ',dist)
                            self.Shift =  0.1*Sphere.r * nplan
                            print( nplan)
                            self.Cylten = Sphere
                            return 1
            return 0 

    def ShiftAll(self):
        self.Shift *= np.random.rand() * 0.1 + 0.9
        if self.shape == 'Cyls':
            self.ShiftAllC()
        if self.shape == 'Sphere':
            self.ShiftAllS()
    def ShiftAllC(self):
        print('All the cylinders are being shifted by Shift ='+str(np.abs(self.Shift)))
        for Cyl in self.Cyls[:]:
            inside1x = 0
            inside1y = 0
            inside1z = 0

            inside2x = 0
            inside2xp = 0
            inside2y = 0
            inside2yp = 0
            inside2z = 0
            inside2zp = 0
            OPn = Cyl[0].OP
            Cyl[0].OP = Cyl[0].OP + np.abs(self.Shift)
            OPnplus1 = Cyl[0].OP
            Cyl[0].p1 = Cyl[0].OP + Cyl[0].e1 * Cyl[0].length/2
            Cyl[0].p2 = Cyl[0].OP - Cyl[0].e1 * Cyl[0].length/2
            
            if OPn[0] < self.sizex-(Cyl[0].L+Cyl[0].marge) and OPn[0] > +Cyl[0].L+Cyl[0].marge:
                inside1x = 1     
            if OPn[1] < self.sizey-(Cyl[0].L+Cyl[0].marge) and OPn[1] > +Cyl[0].L+Cyl[0].marge:
                inside1y = 1     
            if OPn[2] < self.sizez-(Cyl[0].L+Cyl[0].marge) and OPn[2] > +Cyl[0].L+Cyl[0].marge:
                inside1z = 1     
            
            if OPnplus1[0] > self.sizex-(Cyl[0].L+Cyl[0].marge):
                inside2xp = 1     

            if OPnplus1[0] < +Cyl[0].L+Cyl[0].marge:
                inside2x = 1     

            if OPnplus1[1] > self.sizey-(Cyl[0].L+Cyl[0].marge):
                inside2yp = 1     

            if OPnplus1[1] < +Cyl[0].L+Cyl[0].marge:
                inside2y = 1     
            
            if OPnplus1[2] > self.sizez-(Cyl[0].L+Cyl[0].marge):
                inside2zp = 1    
            
            if OPnplus1[2] < +Cyl[0].L+Cyl[0].marge:
                inside2z = 1    

            if  inside1x == 1 and inside2x == 1:
                self.Cyls.append([Cylinder((Cyl[0].OP[0]+self.sizex,Cyl[0].OP[1],Cyl[0].OP[2]),Cyl[0].r, Cyl[0].length, Cyl[0].e1),str(Cyl[1])+'_F'] )   #CylF =
                print(self.Cyls[-1],len(self.Cyls),self.Cyls[-1][1])
            if  inside1x == 1 and inside2xp == 1:
                self.Cyls.append([Cylinder((Cyl[0].OP[0]-self.sizex,Cyl[0].OP[1],Cyl[0].OP[2]),Cyl[0].r, Cyl[0].length, Cyl[0].e1),str(Cyl[1])+'_B'] )   #CylF =
                print(self.Cyls[-1],len(self.Cyls),self.Cyls[-1][1])

            if  inside1y == 1 and inside2y == 1:
                self.Cyls.append([Cylinder((Cyl[0].OP[0],Cyl[0].OP[1]+self.sizey,Cyl[0].OP[2]),Cyl[0].r, Cyl[0].length, Cyl[0].e1),str(Cyl[1])+'_E'] )   #CylE =
                print(self.Cyls[-1],len(self.Cyls),self.Cyls[-1][1])
            if  inside1y == 1 and inside2yp == 1:
                self.Cyls.append([Cylinder((Cyl[0].OP[0],Cyl[0].OP[1]-self.sizey,Cyl[0].OP[2]),Cyl[0].r, Cyl[0].length, Cyl[0].e1),str(Cyl[1])+'_W'] )   #CylE =
                print(self.Cyls[-1],len(self.Cyls),self.Cyls[-1][1])

            if  inside1z == 1 and inside2z == 1:
                self.Cyls.append([Cylinder((Cyl[0].OP[0],Cyl[0].OP[1],Cyl[0].OP[2]+self.sizey),Cyl[0].r, Cyl[0].length, Cyl[0].e1),str(Cyl[1])+'_N'] )   #CylN =
                print(self.Cyls[-1],len(self.Cyls),self.Cyls[-1][1])
            if  inside1z == 1 and inside2zp == 1:
                self.Cyls.append([Cylinder((Cyl[0].OP[0],Cyl[0].OP[1],Cyl[0].OP[2]-self.sizey),Cyl[0].r, Cyl[0].length, Cyl[0].e1),str(Cyl[1])+'_S'] )   #CylN =
                print(self.Cyls[-1],len(self.Cyls),self.Cyls[-1][1])
            
        for Cyl in self.Cyls:
            Cyl[0].InZone = 0
            if Cyl[0].OP[0] < self.sizex+Cyl[0].L+Cyl[0].marge and Cyl[0].OP[0] > -Cyl[0].L-Cyl[0].marge:
                if Cyl[0].OP[1] < self.sizey+Cyl[0].L+Cyl[0].marge and Cyl[0].OP[1] > -Cyl[0].L-Cyl[0].marge:
                    if Cyl[0].OP[2] < self.sizez+Cyl[0].L+Cyl[0].marge and Cyl[0].OP[2] > -Cyl[0].L-Cyl[0].marge:
                        Cyl[0].InZone = 1  

        self.Cyls = [Cyl for Cyl in self.Cyls if Cyl[0].InZone == 1]


    def ShiftAllS(self):
        print('All the Spheres are being shifted by Shift ='+str(np.abs(self.Shift)))
        for Sphere in self.Spheres[:]:
            inside1x = 0
            inside1y = 0
            inside1z = 0

            inside2x = 0
            inside2xp = 0
            inside2y = 0
            inside2yp = 0
            inside2z = 0
            inside2zp = 0
            OPn = Sphere[0].OP
            Sphere[0].OP = Sphere[0].OP + np.abs(self.Shift)
            OPnplus1 = Sphere[0].OP
            
            if OPn[0] < self.sizex-(Sphere[0].r+Sphere[0].marge) and OPn[0] > +Sphere[0].r+Sphere[0].marge:
                inside1x = 1     
            if OPn[1] < self.sizey-(Sphere[0].r+Sphere[0].marge) and OPn[1] > +Sphere[0].r+Sphere[0].marge:
                inside1y = 1     
            if OPn[2] < self.sizez-(Sphere[0].r+Sphere[0].marge) and OPn[2] > +Sphere[0].r+Sphere[0].marge:
                inside1z = 1     
            
            if OPnplus1[0] > self.sizex-(Sphere[0].r+Sphere[0].marge):
                inside2xp = 1     

            if OPnplus1[0] < +Sphere[0].r+Sphere[0].marge:
                inside2x = 1     

            if OPnplus1[1] > self.sizey-(Sphere[0].r+Sphere[0].marge):
                inside2yp = 1     

            if OPnplus1[1] < +Sphere[0].r+Sphere[0].marge:
                inside2y = 1     
            
            if OPnplus1[2] > self.sizez-(Sphere[0].r+Sphere[0].marge):
                inside2zp = 1    
            
            if OPnplus1[2] < +Sphere[0].r+Sphere[0].marge:
                inside2z = 1    

            if  inside1x == 1 and inside2x == 1:
                self.Spheres.append([SphereObj((Sphere[0].OP[0]+self.sizex,Sphere[0].OP[1],Sphere[0].OP[2]),Sphere[0].r),str(Sphere[1])+'_F'] )   #SphereF =
                print(self.Spheres[-1],len(self.Spheres),self.Spheres[-1][1])
            if  inside1x == 1 and inside2xp == 1:
                self.Spheres.append([SphereObj((Sphere[0].OP[0]-self.sizex,Sphere[0].OP[1],Sphere[0].OP[2]),Sphere[0].r),str(Sphere[1])+'_B'] )   #SphereF =
                print(self.Spheres[-1],len(self.Spheres),self.Spheres[-1][1])

            if  inside1y == 1 and inside2y == 1:
                self.Spheres.append([SphereObj((Sphere[0].OP[0],Sphere[0].OP[1]+self.sizey,Sphere[0].OP[2]),Sphere[0].r),str(Sphere[1])+'_E'] )   #SphereE =
                print(self.Spheres[-1],len(self.Spheres),self.Spheres[-1][1])
            if  inside1y == 1 and inside2yp == 1:
                self.Spheres.append([SphereObj((Sphere[0].OP[0],Sphere[0].OP[1]-self.sizey,Sphere[0].OP[2]),Sphere[0].r),str(Sphere[1])+'_W'] )   #SphereE =
                print(self.Spheres[-1],len(self.Spheres),self.Spheres[-1][1])

            if  inside1z == 1 and inside2z == 1:
                self.Spheres.append([SphereObj((Sphere[0].OP[0],Sphere[0].OP[1],Sphere[0].OP[2]+self.sizey),Sphere[0].r),str(Sphere[1])+'_N'] )   #SphereN =
                print(self.Spheres[-1],len(self.Spheres),self.Spheres[-1][1])
            if  inside1z == 1 and inside2zp == 1:
                self.Spheres.append([SphereObj((Sphere[0].OP[0],Sphere[0].OP[1],Sphere[0].OP[2]-self.sizey),Sphere[0].r),str(Sphere[1])+'_S'] )   #SphereN =
                print(self.Spheres[-1],len(self.Spheres),self.Spheres[-1][1])
            
        for Sphere in self.Spheres:
            Sphere[0].InZone = 0
            if Sphere[0].OP[0] < self.sizex+Sphere[0].r+Sphere[0].marge and Sphere[0].OP[0] > -Sphere[0].r-Sphere[0].marge:
                if Sphere[0].OP[1] < self.sizey+Sphere[0].r+Sphere[0].marge and Sphere[0].OP[1] > -Sphere[0].r-Sphere[0].marge:
                    if Sphere[0].OP[2] < self.sizez+Sphere[0].r+Sphere[0].marge and Sphere[0].OP[2] > -Sphere[0].r-Sphere[0].marge:
                        Sphere[0].InZone = 1  

        self.Spheres = [Sphere for Sphere in self.Spheres if Sphere[0].InZone == 1]
            # location des somment millieux 


    def Donne_post_generations(self):
        print('##############################################################################')

        # this is my cylender excluding the ones from periodic copy
        if self.shape == 'Cyls':
            CylsVol = 0
            for Cyl in self.RealCyls:
                CylsVol += math.pi*Cyl[0].r**2*Cyl[0].length
            self.phi  = CylsVol / (self.sizex*self.sizey*self.sizez)
            pp = np.zeros((3,3))
            Ids = [cyl[1].split('_')[0] for cyl in self.Cyls]
            Ids2 = set([cyl[1].split('_')[0] for cyl in self.Cyls])
            for Id in Ids2:
                index = Ids.index(Id)
                
                pp += np.outer(self.Cyls[index][0].e1,self.Cyls[index][0].e1)
            self.OT = pp/len(Ids2)
            print('volume fraction  :', self.phi)
            print( 'Orientation tensor ===>')
            print( '   ')
            print( self.OT)



        if self.shape == 'Sphere':
            SpheresVol = 0
            for Sphere in self.RealSpheres:
                SpheresVol += 4./3. * math.pi*Sphere[0].r**3
            self.phi  = SpheresVol / (self.sizex*self.sizey*self.sizez)
            print ('volume fraction  :', self.phi)
 
    def createCylperio(self,Cyl):
        self.CylsPerios = []
        self.CylsPerios.append([Cylinder((Cyl.OP[0],Cyl.OP[1],Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'N'] )   #CylN =
        self.CylsPerios.append([Cylinder((Cyl.OP[0],Cyl.OP[1],Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'S'] )   #CylS =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1],Cyl.OP[2]),Cyl.r, Cyl.length, Cyl.e1),'F'] )   #CylF =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1],Cyl.OP[2]),Cyl.r, Cyl.length, Cyl.e1),'B'] )   #CylB =
        self.CylsPerios.append([Cylinder((Cyl.OP[0],Cyl.OP[1]+self.sizey,Cyl.OP[2]),Cyl.r, Cyl.length, Cyl.e1),'E'] )   #CylE =
        self.CylsPerios.append([Cylinder((Cyl.OP[0],Cyl.OP[1]-self.sizey,Cyl.OP[2]),Cyl.r, Cyl.length, Cyl.e1),'W'] )   #CylW =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1]+self.sizey,Cyl.OP[2]),Cyl.r, Cyl.length, Cyl.e1),'FE'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1]-self.sizey,Cyl.OP[2]),Cyl.r, Cyl.length, Cyl.e1),'FW'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1]+self.sizey,Cyl.OP[2]),Cyl.r, Cyl.length, Cyl.e1),'BE'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1]-self.sizey,Cyl.OP[2]),Cyl.r, Cyl.length, Cyl.e1),'BW'] )
        self.CylsPerios.append([Cylinder((Cyl.OP[0],Cyl.OP[1]+self.sizey,Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'NW'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0],Cyl.OP[1]-self.sizey,Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'NE'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1],Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'NF'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1],Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'NB'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0],Cyl.OP[1]+self.sizey,Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'SW'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0],Cyl.OP[1]-self.sizey,Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'SE'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1],Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'SF'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1],Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'SB'] ) #Cyl =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1]+self.sizey,Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'NFE'] ) #CylNFE =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1]-self.sizey,Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'NFW'] ) #CylNFW =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1]+self.sizey,Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'NBE'] ) #CylNBE =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1]-self.sizey,Cyl.OP[2]+self.sizez),Cyl.r, Cyl.length, Cyl.e1),'NBW'] ) #CylNBW =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1]+self.sizey,Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'SFE'] ) #CylSFE =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]+self.sizex,Cyl.OP[1]-self.sizey,Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'SFW'] ) #CylSFW =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1]+self.sizey,Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'SBE'] ) #CylSBE =
        self.CylsPerios.append([Cylinder((Cyl.OP[0]-self.sizex,Cyl.OP[1]-self.sizey,Cyl.OP[2]-self.sizez),Cyl.r, Cyl.length, Cyl.e1),'SBW'] ) #CylSBW =
 
    def createSphereperio(self,Sphere):
        self.SpheresPerios = []
        self.SpheresPerios.append([SphereObj((Sphere.OP[0],Sphere.OP[1],Sphere.OP[2]+self.sizez),Sphere.r),'N'] )   #SphereN =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0],Sphere.OP[1],Sphere.OP[2]-self.sizez),Sphere.r),'S'] )   #SphereS =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1],Sphere.OP[2]),Sphere.r),'F'] )   #SphereF =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1],Sphere.OP[2]),Sphere.r),'B'] )   #SphereB =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0],Sphere.OP[1]+self.sizey,Sphere.OP[2]),Sphere.r),'E'] )   #SphereE =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0],Sphere.OP[1]-self.sizey,Sphere.OP[2]),Sphere.r),'W'] )   #SphereW =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1]+self.sizey,Sphere.OP[2]),Sphere.r),'FE'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1]-self.sizey,Sphere.OP[2]),Sphere.r),'FW'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1]+self.sizey,Sphere.OP[2]),Sphere.r),'BE'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1]-self.sizey,Sphere.OP[2]),Sphere.r),'BW'] )
        self.SpheresPerios.append([SphereObj((Sphere.OP[0],Sphere.OP[1]+self.sizey,Sphere.OP[2]+self.sizez),Sphere.r),'NW'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0],Sphere.OP[1]-self.sizey,Sphere.OP[2]+self.sizez),Sphere.r),'NE'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1],Sphere.OP[2]+self.sizez),Sphere.r),'NF'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1],Sphere.OP[2]+self.sizez),Sphere.r),'NB'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0],Sphere.OP[1]+self.sizey,Sphere.OP[2]-self.sizez),Sphere.r),'SW'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0],Sphere.OP[1]-self.sizey,Sphere.OP[2]-self.sizez),Sphere.r),'SE'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1],Sphere.OP[2]-self.sizez),Sphere.r),'SF'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1],Sphere.OP[2]-self.sizez),Sphere.r),'SB'] ) #Sphere =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1]+self.sizey,Sphere.OP[2]+self.sizez),Sphere.r),'NFE'] ) #SphereNFE =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1]-self.sizey,Sphere.OP[2]+self.sizez),Sphere.r),'NFW'] ) #SphereNFW =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1]+self.sizey,Sphere.OP[2]+self.sizez),Sphere.r),'NBE'] ) #SphereNBE =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1]-self.sizey,Sphere.OP[2]+self.sizez),Sphere.r),'NBW'] ) #SphereNBW =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1]+self.sizey,Sphere.OP[2]-self.sizez),Sphere.r),'SFE'] ) #SphereSFE =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]+self.sizex,Sphere.OP[1]-self.sizey,Sphere.OP[2]-self.sizez),Sphere.r),'SFW'] ) #SphereSFW =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1]+self.sizey,Sphere.OP[2]-self.sizez),Sphere.r),'SBE'] ) #SphereSBE =
        self.SpheresPerios.append([SphereObj((Sphere.OP[0]-self.sizex,Sphere.OP[1]-self.sizey,Sphere.OP[2]-self.sizez),Sphere.r),'SBW'] ) #CylSBW =
          
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# import numpy as np
# import mpmath
# import matplotlib
# from mpl_toolkits.mplot3d import Axes3D

# import numpy as np
# from scipy.optimize import curve_fit
# import matplotlib.pyplot as plt
# import scipy.special as sc
# fig = plt.figure()
# ax = Axes3D(fig)

# Cyl1 = Cylinder((2,2,2),0.5,2,(1,0,0))
# size = 100
# Cyl1.pres = 100
# Cyl1.gap = 0.5
# Cyl1.calcul_list(size)
# Cyl2 = Cylinder((2.5,4,2),0.5,2,(0,1,0))
# Cyl2.pres = 100
# Cyl2.gap = 0.5
# Cyl2.calcul_list(size)
# X,Y,Z = [],[],[]

# Cyl1.colisionCheck([[Cyl2,'name']])

# X.extend(Cyl1.X)
# X.extend(Cyl2.X)
# Y.extend(Cyl1.Y)
# Y.extend(Cyl2.Y)
# Z.extend(Cyl1.Z)
# Z.extend(Cyl2.Z)
# Xp1,Yp1,Zp1 = [],[],[]
# Xp2,Yp2,Zp2 = [],[],[]

# Pts1 = []
# Pts1.extend(Cyl1.Pts)
# for pt in Pts1:
#     Xp1.append(pt[0])
#     Yp1.append(pt[1])
# print Xp1
# print Yp1
# Pts1 = []
# Pts1.extend(Cyl1.otherPts)
# for pt in Pts1:
#     Xp2.append(pt[0])
#     Yp2.append(pt[1])
# print Xp2
# print Yp2

# print Cyl1.notok
# print Cyl1.gap
# ax.scatter(Cyl1.X,Cyl1.Y,Cyl1.Z,color='red')
# ax.scatter(Cyl2.X,Cyl2.Y,Cyl2.Z,color='green')
# ax.scatter(Xp1,Yp1,5*np.ones(len(Xp1)),color='red')
# ax.scatter(Xp2,Yp2,5*np.ones(len(Xp2)),color='green')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.show()
