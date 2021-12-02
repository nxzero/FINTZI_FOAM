import AP_class
import math
import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv

class Cylinder():
    def __init__(self,pos,r,length,u):
        self.gap = 0
        self.r = r
        self.OP = np.array(pos)
        self.length = length
        self.e1 = np.array(u)
        self.L = np.sqrt(self.r**2 +  (self.length/2)**2) 
        self.marge = self.L * 0.25
        self.e1 = self.e1/LA.norm(self.e1,2)
        self.AngleChill = np.arctan(2*self.r /self.length)
        ### calce1l de e1 perp 1 ###
        if not abs(self.e1[0]) <= 1e-5:
            self.e2 = np.array([-self.e1[2],0,self.e1[0]])
        elif not abs(self.e1[1]) <= 1e-5:
            self.e2 = np.array([-self.e1[1],self.e1[0],0])
        elif not abs(self.e1[2]) <= 1e-5:
            self.e2 = np.array([0,self.e1[2],self.e1[1]])
        else:
            print("zeros vector !!!!!!!!!!!!!!!!")
        self.e2 = self.e2/LA.norm(self.e2,2)
        self.e3 = np.cross(self.e2 , self.e1)

        self.p1 = self.OP + self.e1 * self.length/2
        self.p2 = self.OP - self.e1 * self.length/2
        self.pres = 1000
        self.InZone = 0
        self.PressureCorrector = 0
        self.ppRef =  0

        self.gamma = list(np.linspace(0,math.pi*2,self.pres))
        self.i = 0
        self.gapbetweenplanes = self.r / 30.

    def calcul_planKey(self,sizex,sizey,sizez):
        self.sizex = sizex
        self.sizey = sizey
        self.sizez = sizez
        # i should only calculs points if i keep the cylinder
        
        self.list_of_point = []
        self.pts = {}
        self.keys =[
            'H','B',
            'Ixb','Ixt','Iyb','Iyt','Izb','Izt',
            'HIxb','HIxt','HIyb','HIyt','HIzb','HIzt',
            'BIxb','BIxt','BIyb','BIyt','BIzb','BIzt'
            ]
        for key in self.keys:
            self.pts[key] = {
                'X':[],
                'Y':[],
                'Z':[]
            }
        # definition des plans
        self.nxb = np.array([-1,0,0])
        self.nxt = np.array([1,0,0])
        self.nyb = np.array([0,-1,0])
        self.nyt = np.array([0,1,0])
        self.nzb = np.array([0,0,-1])
        self.nzt = np.array([0,0,1])

        self.origb = np.array([0,0,0])
        self.origt = np.array([self.sizex,self.sizey,self.sizez])

        self.plans = {
            'Ixb':[self.nxb,self.origb],
            'Ixt':[self.nxt,self.origt],
            'Iyb':[self.nyb,self.origb],
            'Iyt':[self.nyt,self.origt],
            'Izb':[self.nzb,self.origb],
            'Izt':[self.nzt,self.origt]
        }

        self.plansKeys  = set()

        for key,l in [['H',-self.length],['B',self.length]]:
            self.PC = self.e1 * l/2
            self.OC = self.OP + self.PC
            for g in self.gamma:
                self.CM = self.r * (math.cos(g)*self.e2 + math.sin(g) * self.e3)
                self.OM = self.OC + self.CM
                self.list_of_point.append(self.OM)
                if self.OM[0] < self.sizex and self.OM[0] > 0:
                    if self.OM[1] < self.sizey and self.OM[1] > 0:
                        if self.OM[2] < self.sizez and self.OM[2] > 0:
                            self.pts[key]['X'].append(self.OM[0])
                            self.pts[key]['Y'].append(self.OM[1])
                            self.pts[key]['Z'].append(self.OM[2])
                if self.OM[0] > self.sizex:
                    self.plansKeys.add('Ixt')
                if self.OM[0] < 0:
                    self.plansKeys.add('Ixb')
                if self.OM[1] > self.sizey:
                    self.plansKeys.add('Iyt')
                if self.OM[1] < 0:
                    self.plansKeys.add('Iyb')
                if self.OM[2] > self.sizez:
                    self.plansKeys.add('Izt')
                if self.OM[2] < 0:
                    self.plansKeys.add('Izb') 

                

        if (self.pts['H']['X'] == [] and self.pts['B']['X'] == [] and len(self.plansKeys) == 1):
            self.plansKeys = []


        pKeytmp = list(self.plansKeys)[:]
        if len(self.pts['H']['X']) != 0 and len(self.pts['H']['X']) != self.pres:
            for key2 in pKeytmp:
                if key2[1] == 'y' and abs(min(self.pts['H']['Y'])) < self.gapbetweenplanes:
                    self.plansKeys.add('H'+key2)
                if key2[1] == 'y' and abs(max(self.pts['H']['Y']) - self.sizey) < self.gapbetweenplanes:
                    self.plansKeys.add('H'+key2)
                if key2[1] == 'x' and abs(min(self.pts['H']['X'])) < self.gapbetweenplanes:
                    self.plansKeys.add('H'+key2)
                if key2[1] == 'x' and abs(max(self.pts['H']['X']) - self.sizex) < self.gapbetweenplanes:
                    self.plansKeys.add('H'+key2)
                if key2[1] == 'z' and abs(min(self.pts['H']['Z'])) < self.gapbetweenplanes:
                    self.plansKeys.add('H'+key2)
                if key2[1] == 'z' and abs(max(self.pts['H']['Z']) - self.sizez) < self.gapbetweenplanes:
                    self.plansKeys.add('H'+key2)


        if len(self.pts['B']['X']) != 0 and len(self.pts['B']['X']) != self.pres:
            for key2 in [keys for keys in pKeytmp if keys[0] != 'H' or keys[0] != 'B']:
                if key2[1] == 'y' and abs(min(self.pts['B']['Y'])) < self.gapbetweenplanes:
                    self.plansKeys.add('B'+key2)
                if key2[1] == 'y' and abs(max(self.pts['B']['Y']) - self.sizey) < self.gapbetweenplanes:
                    self.plansKeys.add('B'+key2)
                if key2[1] == 'x' and abs(min(self.pts['B']['X'])) < self.gapbetweenplanes:
                    self.plansKeys.add('B'+key2)
                if key2[1] == 'x' and abs(max(self.pts['B']['X']) - self.sizex) < self.gapbetweenplanes:
                    self.plansKeys.add('B'+key2)
                if key2[1] == 'z' and abs(min(self.pts['B']['Z'])) < self.gapbetweenplanes:
                    self.plansKeys.add('B'+key2)
                if key2[1] == 'z' and abs(max(self.pts['B']['Z']) - self.sizez) < self.gapbetweenplanes:
                    self.plansKeys.add('B'+key2)

    def calcul_list(self):
        self.OI = 'none'
        for key in [keys for keys in self.plansKeys if keys[0]=='I']:
            
            #normal distance to the plan
            n = self.plans[key][0]
            Ip = self.plans[key][1]
            
            #(pts on plan - OP) dot the normal
            dist = np.dot((Ip - self.OP),n)
            if dist < self.length:
                if  np.dot(self.e1,n) != 0:
                    t = dist / np.dot(self.e1,n)
                    self.OI = self.OP + self.e1*t
                else:
                    t = dist
                    self.OI = self.OP + n*t
                #this is the point of intersection
                #let's define the basis in the plane
                ep2 = np.cross(n,self.e1)
                if LA.norm(ep2,2) != 0:
                    ep2 = ep2 / LA.norm(ep2,2)
                else:
                    ep2 = np.cross(n,self.e2)
                    ep2 = ep2 / LA.norm(ep2,2)
                ep1 = np.cross(ep2,n)
                ep1 = ep1 / LA.norm(ep1,2)
                #find theta
                # if np.dot(ep1,self.e1) != 1 or np.dot(ep1,self.e1) != -1:
                if  np.abs(1-np.dot(ep1,self.e1)) >= 1e-10:
                    theta = np.arccos(np.dot(ep1,self.e1))
                else:
                    theta = 0

                # else:
                #     theta = 0
                if np.sin(theta) != 0:
                    h = self.r / np.sin(theta)
                else:
                    h = self.length * 50 #look like straight lines
                if theta <= 0.1:
                    self.gamma = list(np.linspace(0,math.pi*2,self.pres*50))
                if theta == 0:
                    self.gamma = list(np.linspace(0,math.pi*2,self.pres*200))
                
                self.epsilon = self.sizex*1e-6

                for g in self.gamma:
                    # repere local
                    x = h * np.cos(g)
                    y = self.r* np.sin(g)
                    # repere global
                    ptss = self.OI + x*ep1 +y*ep2 - n * self.epsilon*2
                    if (LA.norm(ptss-self.OP,2)) <= self.L:
                        if ptss[0] < self.sizex and ptss[0] > 0:
                            if ptss[1] < self.sizey and ptss[1] > 0:
                                 if ptss[2] < self.sizez and ptss[2] > 0:
                                    self.pts[key]['X'].append(ptss[0])
                                    self.pts[key]['Y'].append(ptss[1])
                                    self.pts[key]['Z'].append(ptss[2])
                self.gamma = list(np.linspace(0,math.pi*2,self.pres))

        #intersection plan-plan        






        ##########""

            
        line = list(np.linspace(-1,1,self.pres))

        ###### pts d'intersection avec les bords ############
        # 6 bors to check for each two cicrles of the cylinder
        # easy
        for key in [keys for keys in self.plansKeys if keys[0:2]=='HI' or keys[0:2]=='BI']:

            #normal distance to the plan
            n = self.plans[key[1:]][0]
            I = self.plans[key[1:]][1]
            if key[0] == 'B':
                l = self.length
            if key[0] == 'H':
                l = -self.length
            C = self.OP + self.e1 * l/2
            if LA.norm(np.cross(self.e1,n),2) != 0:
                n2 = np.cross(n,self.e1)/LA.norm(np.cross(self.e1,n),2)
                n3 = np.cross(n2,n)
            else:
                continue
            #e1 in the local frame (n n2 n3)
            el = np.array([ np.dot(self.e1,n3) , np.dot(self.e1,n) ])
            el = el/LA.norm(el,2)
            #the normal of the plan in the local frame
            nl = np.array([0,1])


            #C in the local frame
            Cl = np.array([ np.dot(C,n3), np.dot(C,n) ])
            Il = np.array([ np.dot(I,n3), np.dot(I,n) ])


            Bmat = np.array([np.dot(Cl,el),np.dot(Il,nl)])
            Amat = np.array([
                [el[0],el[1]],
                [nl[0],nl[1]]
            ])


            Xl = np.dot(inv(Amat),Bmat)

            X = Xl[0]*n3 + Xl[1]*n + n2*np.dot(C,n2)
   
            for g in line:
                ptss = X+n2*g*self.r - n * self.epsilon
                if (LA.norm(ptss-self.OP,2)) <= self.L:
                    if ptss[0] < self.sizex and ptss[0] > 0:
                        if ptss[1] < self.sizey and ptss[1] > 0:
                            if ptss[2] < self.sizez and ptss[2] > 0:
                                self.pts[key]['X'].append(ptss[0])
                                self.pts[key]['Y'].append(ptss[1])
                                self.pts[key]['Z'].append(ptss[2])

            
        self.M = {}
        X,Y,Z = list(),list(),list()
        for key in self.keys:
            self.M[key] = list(zip(self.pts[key]['X'],self.pts[key]['Y'],self.pts[key]['Z']))

    def rearrangeList(self):
        for key in self.keys:
            index = []
            if len(list(self.M[key])) > 3:
                dist1  =  np.array(self.M[key][0]) -np.array(self.M[key][1]) 
                dist2  =  np.array(self.M[key][1]) -np.array(self.M[key][2]) 
                dist = min(LA.norm(dist1,2),LA.norm(dist2,2))
                for i in range(len(self.M[key])):
                    D = np.array(self.M[key][i]) -np.array(self.M[key][(i+1)%len(self.M[key])])
                    if LA.norm(D,2) > 3*dist:
                        index.append(i)
        
            if len(index) == 1:
                i = index[0]
                self.M[key] = self.M[key][(i+1)%len(self.M[key]):] + self.M[key][:(i+1)%len(self.M[key])]

            if len(index) == 2:
                i = index[0]
                keybis = key +'_bis'
                self.keys.append(keybis)
                self.M[keybis] = self.M[key][(i+1)%len(self.M[key]):]
                self.M[key] = self.M[key][:(i+1)%len(self.M[key])]
            if len(index) == 3:
                i = index[0]
                i2 = index[1]
                keybis = key +'_bis'
                keyter = key +'_ter'
                self.keys.append(keybis)
                self.keys.append(keyter)
                self.M[keyter] = self.M[key][(i2+1)%len(self.M[key]):]
                self.M[keybis] = self.M[key][(i+1)%len(self.M[key]):(i2+1)%len(self.M[key])]
                self.M[key] = self.M[key][:(i+1)%len(self.M[key])]

    def colisionCheck(self,others):
        self.notok = 0
        self.mmplan = 0
        others = [other[0] for other in others]
        for other in others:
            if not self == other:
                if LA.norm(self.OP - other.OP,2) >= (self.L + other.L +self.gap):
                    self.notok = 0
                    continue
                ### rajouter le critere shpere sphere
                self.normal = np.cross(self.e1,other.e1) #normal au plan commun
                if LA.norm(self.normal,2) != 0:
                    self.normal = self.normal / LA.norm(self.normal,2)
                else:
                    vec = other.OP - self.OP
                    self.normal = np.cross(self.e1,vec)
                    if LA.norm(self.normal,2) != 0:
                        self.normal = self.normal / LA.norm(self.normal,2)
                

                self.dvec = self.OP - other.OP
                self.d = LA.norm(self.dvec,2)
                self.dmin = np.dot(self.dvec,self.normal)
                if LA.norm(self.OP - other.OP,2) <= (self.r + other.r +self.gap):
                    self.notok = 1
                    break

                if abs(self.dmin) <= (self.r + other.r +self.gap):
                    self.mmplan = 1
                    #define a local planar basis :
                    self.e_perp = np.cross(self.e1,self.normal)
                    #define verticis in the 2D geometry
                    self.Pts = list(range(4))
                    self.Pts[0] =np.array([  self.r+self.gap/2., - (self.length/2+self.gap/2.)])
                    self.Pts[1] =np.array([  self.r+self.gap/2.,   self.length/2+self.gap/2.])
                    self.Pts[2] =np.array([- (self.r+self.gap/2.), +  self.length/2+self.gap/2.])
                    self.Pts[3] =np.array([- (self.r+self.gap/2.), -  (self.length/2+self.gap/2.)])
                    #define segment
                    self.Seg = list(range(4))
                    self.Seg[0] = self.Pts[1] - self.Pts[0]
                    self.Seg[1] = self.Pts[2] - self.Pts[1]
                    self.Seg[2] = self.Pts[3] - self.Pts[2]
                    self.Seg[3] = self.Pts[0] - self.Pts[3]
                    # deifne verticis for the other
                    other.e_perp = np.cross(self.normal,other.e1)
                    #define verticis
                    other.br_GB = (other.r+self.gap) * other.e_perp - (other.length/2+self.gap) * other.e1
                    other.tr_GB = (other.r+self.gap) * other.e_perp + (other.length/2+self.gap) * other.e1
                    other.tl_GB = - (other.r+self.gap) * other.e_perp + (other.length/2+self.gap) * other.e1
                    other.bl_GB = - (other.r+self.gap) * other.e_perp - (other.length/2+self.gap) * other.e1

                    #Project on the local basis
                    other.Pts = list(range(4))
                    other.Pts[0] =np.array([np.dot(other.br_GB,self.e_perp), np.dot(other.br_GB,self.e1)])
                    other.Pts[1] =np.array([np.dot(other.tr_GB,self.e_perp), np.dot(other.tr_GB,self.e1)])
                    other.Pts[2] =np.array([np.dot(other.tl_GB,self.e_perp), np.dot(other.tl_GB,self.e1)])
                    other.Pts[3] =np.array([np.dot(other.bl_GB,self.e_perp), np.dot(other.bl_GB,self.e1)])
                    # add the offset :
                    self.OPplan = other.OP - self.OP
                    self.OPplan = np.array([np.dot(self.OPplan,self.e_perp), np.dot(self.OPplan,self.e1)])
                    self.crit = self.L + other.L
                    if LA.norm(self.OPplan,2) > (self.crit+self.gap):
                        continue
                    # projection sur le plan
                    other.Pts[0] = other.Pts[0] + self.OPplan
                    other.Pts[1] = other.Pts[1] + self.OPplan
                    other.Pts[2] = other.Pts[2] + self.OPplan
                    other.Pts[3] = other.Pts[3] + self.OPplan

                    self.otherPts = other.Pts

                    other.Seg = list(range(4))
                    other.Seg[0] = other.Pts[1] - other.Pts[0]
                    other.Seg[1] = other.Pts[2] - other.Pts[1]
                    other.Seg[2] = other.Pts[3] - other.Pts[2]
                    other.Seg[3] = other.Pts[0] - other.Pts[3]
                    # Verification de l intersection des droites
                    for i in list(range(4)):
                        A = self.Pts[i]
                        B = self.Pts[(i+1)%4]
                        AB = B - A
                        for j in list(range(4)):
                            C = other.Pts[j]
                            D = other.Pts[(j+1)%4]
                            CD = D - C
                            M = np.array([[AB[0],-CD[0]],[AB[1],-CD[1]]])
                            E = C - A
                            try:
                                invM = inv(M)
                                X = np.dot(invM,E)
                            except:
                                X = [-1,-1]
                            if X[0] < 1 and X[0] > 0 and X[1] < 1 and X[1] > 0:
                                self.notok = 1
                                break
                        else:
                            continue
                        break
                    else:
                        continue
                    break

