import numpy as np
import numpy.linalg as LA
import pandas as pd
import math
class Bubble:
    """The goal is to gather each path of this bubble"""
    def __init__(self,tag: int):
        self.tag = tag
        self.D = 1.             # Define the diametre with respect to the volume since bubbles can merges 
        self.rmax = 3.
        self.dr = self.D/5.
        self.Rs = np.arange(self.D/100,self.rmax,self.dr)
        self.x = list()
        self.y = list()
        self.vx = list()
        self.vy = list()
        self.vol = list()
        self.t = list()
        self.i = list()
        self.g_r = list()
        self.contact = list()
        self.dist_tag = list()  # list of the distance between one this drop and the others referanced by the tag
        
    def assigne_0(self):
        i = int(0)
        t = float(0)
        x = float(0)
        y = float(0)
        vx = float(0)
        vy = float(0)
        vol = float(0)
        self.x.append(x)
        self.y.append(y)
        self.vx.append(vx)
        self.vy.append(vy)
        self.vol.append(vol)
        self.t.append(t)
        self.i.append(i)
        self.contact(0)
        self.g_r.append(np.zeros(len(self.Rs)))    
        self.dist_tag.append(np.zeros(len(self.Rs)))
        
    def assigne(self,bstate: pd.DataFrame):
        i = int(bstate.i.values)
        t = float(bstate.t.values)
        x = float(bstate.x.values)
        y = float(bstate.y.values)
        vx = float(bstate.vx.values)
        vy = float(bstate.vy.values)
        vol = float(bstate.vol.values)
        self.x.append(x)
        self.y.append(y)
        self.vx.append(vx)
        self.vy.append(vy)
        self.vol.append(vol)
        self.t.append(t)
        self.i.append(i)
        self.g_r.append(list())    
        self.dist_tag.append(list())
        self.contact.append(0)
        
    def assigne_list(self,Pos: pd.DataFrame):
        self.x = Pos.x[Pos.tag.eq(self.tag)].values
        self.y = Pos.y[Pos.tag.eq(self.tag)].values
        self.vx = Pos.vx[Pos.tag.eq(self.tag)].values
        self.vy = Pos.vy[Pos.tag.eq(self.tag)].values
        self.vol = Pos.vol[Pos.tag.eq(self.tag)].values
        self.t = Pos.t[Pos.tag.eq(self.tag)].values
        self.i = Pos.i[Pos.tag.eq(self.tag)].values
        

    

    def dist_perio(self,other,step: int,size: float):
        pos_b = np.array([  self.x[step]  ,  self.y[step]  ])
        other_b = np.array([  other.x[step]  ,  other.y[step]  ])
        if(pos_b[0]>0):
            x_centered = (other_b[0]<pos_b[0]-size/2)*size+np.array(other_b[0])
        else:
            x_centered = (other_b[0]>pos_b[0]+size/2)*(-size)+np.array(other_b[0])
        if(pos_b[1]>0):
            y_centered = (other_b[1]<pos_b[1]-size/2)*size+np.array(other_b[1])
        else:
            y_centered = (other_b[1]>pos_b[1]+size/2)*(-size)+np.array(other_b[1])
            
        pos_centered = np.array([x_centered,y_centered])
        #Now we can calculate the diff between our pts and the centered others
        return LA.norm(pos_centered - pos_b,2)

    
    def mean_time_values(self):
        self.mean_x = sum(self.x)/len(self.x)
        self.mean_y = sum(self.y)/len(self.y)
        self.mean_vx = sum(self.vx)/len(self.vx)
        self.mean_vy = sum(self.vy)/len(self.vy)
        self.mean_vol = sum(self.vol)/len(self.vol)
        self.mean_g_r = sum(np.array(self.g_r))/len(self.g_r)
        
    
    def div(self,other: int):
        self.x = [x/other for x in self.x]
        self.y = [x/other for x in self.y]
        self.vx = [x/other for x in self.vx]
        self.vy = [x/other for x in self.vy]
        self.vol = [x/other for x in self.vol]
        self.g_r = [list(np.array(x)/other) for x in self.g_r]
        return self
        
    def __add__(self,other):
        self.x = [x+y for y,x in list(zip(other.x,self.x))]
        self.y = [x+y for y,x in list(zip(other.y,self.y))]
        self.vx = [x+y for y,x in list(zip(other.vx,self.vx))]
        self.vy = [x+y for y,x in list(zip(other.vy,self.vy))]
        self.vol = [x+y for y,x in list(zip(other.vol,self.vol))]
        self.g_r = [list(np.array(x)+np.array(y)) for y,x in list(zip(other.g_r,self.g_r))]
        return self
        
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)