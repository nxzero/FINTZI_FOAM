import AP_class
import math
import numpy as np
import scipy as sp
from numpy import linalg as LA
from cylinder import *

class SphereObj(Cylinder):
    def __init__(self, pos,r):
        Cylinder.__init__(self,pos,r,r,[1,0,0])
        self.marge = self.r * 0.25

    def calcul_planKey(self,sizex,sizey,sizez):
        self.sizex = sizex
        self.sizey = sizey
        self.sizez = sizez
        # i should only calculs points if i keep the cylinder
        gamma = list(np.linspace(0,math.pi*2,self.pres))
        self.list_of_point = []
        self.pts = {}
        keys =['Ixb','Ixt','Iyb','Iyt','Izb','Izt']
        for key in keys:
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
            'Ixb':[self.nxb,self.origb,self.nyb,self.nzb],
            'Ixt':[self.nxt,self.origt,self.nyt,self.nzt],
            'Iyb':[self.nyb,self.origb,self.nzb,self.nxb],
            'Iyt':[self.nyt,self.origt,self.nzt,self.nxt],
            'Izb':[self.nzb,self.origb,self.nxb,self.nyb],
            'Izt':[self.nzt,self.origt,self.nxt,self.nyt]
        }

        self.plansKeys  = set()


        if np.abs(np.dot(self.OP,self.plans['Ixt'][0])) + self.r >= self.sizex:
            self.plansKeys.add('Ixt')
        if np.abs(np.dot(self.OP,self.plans['Ixb'][0])) - self.r <= 0:
            self.plansKeys.add('Ixb')
        if np.abs(np.dot(self.OP,self.plans['Iyt'][0])) + self.r >= self.sizey:
            self.plansKeys.add('Iyt')
        if np.abs(np.dot(self.OP,self.plans['Iyb'][0])) - self.r <= 0:
            self.plansKeys.add('Iyb')
        if np.abs(np.dot(self.OP,self.plans['Izt'][0])) + self.r >= self.sizez:
            self.plansKeys.add('Izt')
        if np.abs(np.dot(self.OP,self.plans['Izb'][0])) - self.r <= 0:
            self.plansKeys.add('Izb')  

    def calcul_list(self):
        for key in self.plansKeys:

            #normal distance to the plan
            n = self.plans[key][0]
            n1 = self.plans[key][2]
            n2 = self.plans[key][3]
            Ip = self.plans[key][1]
            
            #(pts on plan - OP) dot the normal
            dist = np.dot((Ip - self.OP),n)
            petit_r = np.sqrt(np.abs(self.r**2-dist**2))
            for g in gamma:
                # repere local
                x = petit_r * np.cos(g)
                y = petit_r * np.sin(g)
                # repere global
                self.epsilon = self.sizex*1e-6
                ptss = self.OP+x*n1 +y*n2 + n * (dist-self.epsilon)
                if ptss[0] < self.sizex and ptss[0] > 0:
                    if ptss[1] < self.sizey and ptss[1] > 0:
                        if ptss[2] < self.sizez and ptss[2] > 0:
                            self.pts[key]['X'].append(ptss[0])
                            self.pts[key]['Y'].append(ptss[1])
                            self.pts[key]['Z'].append(ptss[2])
        
        self.M = {}
        X,Y,Z = list(),list(),list()
        for key in keys:
            self.M[key] = zip(self.pts[key]['X'],self.pts[key]['Y'],self.pts[key]['Z'])

    def rearrangeList(self):
        Cylinder.rearrangeList(self)


    def colisionCheck(self,others):
        from numpy.linalg import inv
        self.notok = 0
        self.mmplan = 0
        others = [other[0] for other in others]
        for other in others:
            if not self == other:
                if LA.norm(self.OP - other.OP,2) <= (self.r + other.r +self.gap):
                    self.notok = 1
                    break

