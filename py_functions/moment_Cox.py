
import numpy as np
import scipy.special as sc
import math
# Couple proposee par Khayat  Cox
gamma = 0.5772156649

def T_simple(Re,Chi,Theta):
    import math
    Theta = np.array(Theta) /360 * 2*math.pi
    T = 5.*math.pi/12. *Chi*Re / np.log(Chi)**2 * np.sin(2*Theta)
    return T

def T_ad(Re,Chi,Theta):
    import math
    Theta = np.array(Theta) /360 * 2*math.pi
    Re = np.array(Re)
    Chi = np.array(Chi)
    A = 5.*math.pi/12. * Re/(1.+Chi*Re**1.1)**0.5
    B = Chi/np.log(Chi)**2
    C = (13.5-30.*Re**(1./2.))/(Chi*np.log(Chi)**3)*np.exp(-0.7*Re)
    return A * (B+C) * np.sin(2*Theta)

def Zg(x):
    x = np.array(x)
    Z = 2./x * (1.+(np.exp(-x)-1.)/x)
    return Z
def Cg(X):
    return (sc.expn(1,X)+np.log(X)+gamma)/X

def F_g(ReL,Theta):
    import math
    Theta = np.array(Theta)
    Theta = Theta /360. * 2.*math.pi
    X =ReL*(1.-np.cos(Theta))
    Y =ReL*(1.+np.cos(Theta))
    B = Zg(X) -Cg(X) + Zg(Y) -Cg(Y)
    return (np.cos(Theta) * B + Zg(Y) - Zg(X)) * np.sin(Theta)

def T_CC(Re,Chi,Theta) :
    import math
    Re = np.array(Re)
    Chi = np.array(Chi)
    ReL = Re * Chi / 2.
    A = - 2.* math.pi/np.log(Chi)**2.
    return  A * F_g(ReL,Theta)

def T_CC_myfit(Re,Chi,Theta):#
    try:
        if Re == 0:
            return 0
        if Theta == 0:
            Theta = 0.001
        if Theta == 90:
            Theta = 89.99
    except:
        1
    C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12 = 0.64322646,   0.34011955,  -0.12984728,  30.        ,-0.7       ,   0.5,   0., 0. , 0.,   3.,   1.,   1.5  
    import math
    Theta = np.array(Theta) /360 * 2*math.pi
    Re = np.array(Re)
    Chi = np.array(Chi)
    R = C1*Re**C2+C3
    R2 = (Chi - C4)**2*C5 * np.tanh(Re*C10-C6)
    A = C7 /(np.log(Chi)**1)
    A2 = C8 /(Chi*np.log(Chi)**2)
    A3 = C9 /(Chi*np.log(Chi)**3)
    A4 = C11 /(Chi*np.log(Chi)**4)
    B = C8/(Chi*Re**2*C9*np.log(Chi)**4)
    C = C11*np.exp(-Re*C12)*Chi

    # B =  C8*Chi/np.log(Chi)
    return  R*np.sin(2*Theta)

### Forces Cox ###########################################

def F_90_const(X,C1,C2,C3,C4,C6,C7,C9,C10,C11,C12):
    import math
    Re,Chi = X
    f_perp =sc.expn(1,Re*Chi/2)+np.log(Re*Chi/2)-(np.exp(-Re*Chi/2)-1.)/(Re*Chi/2)+gamma-1./2-np.log(4)
    my_coef  = (1-np.exp(C1*(Chi-C2)))
    my_second_coef =C3 *np.exp(-C4/Chi)  *Re**C12
    my_third_coef = C6 *np.exp(-C7/Chi) *Re**C11
    F_perp2 = 4./3*(Chi)/Chi**(1./3)*(2./3)**(1./3)*(1./np.log(2*Chi) + (1./2-np.log(2)+(f_perp-1./2+2*np.log(2))*my_coef)/np.log(2*Chi)**2 +  (1.- np.exp(-0.01*Chi**2.5))*( my_second_coef *0.21484/np.log(2*Chi)**3 +  (my_third_coef * 0.38735)/np.log(2*Chi)**4))
    return F_perp2

def F_90(Re,Chi):
    p0 = np.array([ -1.22840432e+01,   1.93973798e+00,   3.49519674e+01,2.44794412e+01,   1.01967917e+00,  -3.29851411e-02,3.00000000e+00,   1.00000000e+00,   1.52080416e+00,1.03427954e+00])
    return F_90_const((Re,Chi),*p0)

def F_0_KC(Re,Chi):
    import math
    Re = np.array(Re)
    Chi = np.array(Chi)
    f_A = (sc.expn(1,Chi*Re) + np.log(Chi*Re) - np.exp(-Chi*Re) + 1 +gamma)/(Chi*Re)
    f_B = sc.expn(1,Chi*Re) + np.log(Chi*Re) + gamma -2.
    f_para = 1./2. * (f_A + f_B)
    A = 1./np.log(Chi)
    B = (0.80685 + f_para -np.log(2.))* 1./np.log(Chi)**2
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B)-1.4/(Chi-0.5)**1.75
    return F_para2
    
def F_0(Re,Chi):
    import math
    Re = np.array(Re)
    Chi = np.array(Chi)
    f_A = (sc.expn(1,Chi*Re) + np.log(Chi*Re) - np.exp(-Chi*Re) + 1 +gamma)/(Chi*Re)
    f_B = sc.expn(1,Chi*Re) + np.log(Chi*Re) + gamma -2.
    f_para = 1./2. * (f_A + f_B)
    f_3 = (Chi*Re)**(0.07*Chi**0.5)
    f_4 = (Chi* Re)**(0.03*Chi**0.9)
    A = 1./np.log(2*Chi)
    B = (0.80685 + f_para )* 1./np.log(2*Chi)**2
    C = (0.82854 + f_para  *f_3 )* 1./np.log(2*Chi)**3
    D = (1.45243  + f_para * f_4)* 1./np.log(2*Chi)**4
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D)-1.4/(Chi-0.5)**1.75
    return F_para2

def F_0_myfit(Re,Chi):
    C1,C2,C3,C4,C5,C6 = 0.23906317, -0.02625534,  0.35945921, -0.17541541,  2.17545436,-1.58283729
    Re = np.array(Re)
    Chi = np.array(Chi)
    f_A = (sc.expn(1,Chi*Re) + np.log(Chi*Re) - np.exp(-Chi*Re) + 1 +gamma)/(Chi*Re)
    f_B = sc.expn(1,Chi*Re) + np.log(Chi*Re) + gamma -2.
    f_para = 1./2. * (f_A + f_B)
    f_3 = C5*(Chi*Re)**(C1*Chi**C2)
    f_4 = C6*(Chi* Re)**(C3*Chi**C4)
    A = 1./np.log(2*Chi)
    B = (0.80685 + f_para )* 1./np.log(2*Chi)**2
    C = (0.82854 + f_para  *f_3 )* 1./np.log(2*Chi)**3
    D = (1.45243  + f_para * f_4)* 1./np.log(2*Chi)**4
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D)-1.4/(Chi-0.5)**1.75
    return F_para2


def F_0_modi(Chi):
    import math
    Chi = np.array(Chi)
    f_3 = 1#(Chi*Re)**(0.07*Chi**0.5)
    f_4 = 1#(Chi* Re)**(0.03*Chi**0.9)
    A = 1./np.log(2*Chi)
    B = 0.80685* 1./np.log(2*Chi)**2
    C = 0.82854* 1./np.log(2*Chi)**3
    D = 1.45243* 1./np.log(2*Chi)**4
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D)-1.4/(Chi-0.5)**1.75
    return F_para2

def F_04th(Chi):
    import math
    Chi = np.array(Chi)
    a_1 = 1
    a_2 = 0.80685
    a_3 = 0.82854
    a_4 = 1.45243
    A = 1./np.log(2*Chi)
    B = 0.80685* 1./np.log(2*Chi)**2
    C = 0.82854* 1./np.log(2*Chi)**3
    D = 1.45243* 1./np.log(2*Chi)**4
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D)
    return F_para2
############## para fir JL ##################
def F_para(Re,Chi,theta):
    theta = np.array(theta) /360 * 2*math.pi
    return F_0(Re,Chi) * np.cos(theta)

def F_para_modi(Re,Chi,theta):
    theta = np.array(theta) /360 * 2*math.pi
    return (F_0(Re,Chi) - 1./16. *np.sin(theta)**2  *Chi*Re / (np.log(Chi)**2))* np.cos(theta)

############ para my fits #####################

def F_para_myfit(Re,Chi,theta):
    theta = np.array(theta) /360 * 2*math.pi
    return F_0_myfit(Re,Chi) * np.cos(theta)

def F_para_modi_myfit(Re,Chi,theta):
    try:    
        if Re == 0:
            return 0
        if theta == 0:
            theta = 0.0001
        if theta == 90:
            theta = 89.999
    except:
        1
    C1, C2 = 0.00820158,  1.48790879
    theta = np.array(theta) /360 * 2*math.pi
    return (F_0_myfit(Re,Chi) + C1 *np.sin(theta*C2)**2  *Re*Chi )* np.cos(theta)


###########perp myfits ####################


def F_perp(Re,Chi,theta):
    theta = np.array(theta) /360 * 2*math.pi
    return F_90(Re,Chi) * np.sin(theta)

def F_perp_Modi(Re,Chi,theta):
    theta = np.array(theta) /360 * 2*math.pi
    return (F_90(Re,Chi) + 1./16. *np.cos(theta)**2  *Chi*Re / (np.log(Chi)**2)) * np.sin(theta)

def F_perp_modi_myfit(Re,Chi,theta ):
    try:    
        if Re == 0:
            return 0
        if theta == 0:
            theta = 0.0001
        if theta == 90:
            theta = 89.999
    except:
        1
    C1,C2 = 0.01297534,  1.11691914
    theta = np.array(theta) /360 * 2*math.pi
    return (F_90(Re,Chi) - C1 *np.cos(theta*C2)**2  *Chi*Re) * np.sin(theta)

###########" avec les valeurs de la simu ###########
def F_paraD(Re,Chi,theta,F0):
    theta = np.array(theta) /360. * 2.*math.pi
    return F0 * np.cos(theta)

def F_perpD(Re,Chi,theta,F90):
    theta = np.array(theta) /360. * 2.*math.pi
    return F90 * np.sin(theta)

def F_perpD_Modi(Re,Chi,theta,F90):
    theta = np.array(theta) /360. * 2.*math.pi
    return (F90 - 1./16. *np.cos(theta)**2  *Chi*Re / (np.log(Chi)**2)) * np.sin(theta)

def F_paraD_modi(Re,Chi,theta,F0):
    theta = np.array(theta) /360. * 2.*math.pi
    return (F0 + 1./16. *np.sin(theta)**2  *Chi*Re / (np.log(Chi)**2))* np.cos(theta)

####################################################
########### Forces Cox 2nd ordre ##################

def C(x):
    return sc.expn(1,x)+np.log(x) + gamma - x
def Z(x):
    return (np.exp(-x) - 1)/x
def B(x):
    return 1./2. * (sc.expn(1,x) - Z(x))

def F_D(ReL,Theta):#theta en radians
    X = ReL*(1.-np.cos(Theta))
    Y = ReL*(1.+np.cos(Theta))
    CA = np.cos(Theta)**2 / (2.*ReL) * (C(X) + C(Y)) + 3.*np.cos(Theta)**2 - 2.
    CB = B(X) + B(Y) + 1./2.*np.log(1.-np.cos(Theta)**2) + gamma + np.log(ReL / 4.)
    return CA*1./(2.*(2.-np.cos(Theta)**2)) + CB


def D(Re,Chi,Theta): # diverge en 0
    Re = np.array(Re)
    Chi = np.array(Chi)
    ReL = Re* Chi/2.
    try:
        if max(Theta) >= 3: 
            Theta = Theta / 360. *2. *math.pi
    except:
        if Theta >= 3:
            Theta = Theta / 360. *2. *math.pi

    A = - 4. * math.pi * (2. - np.cos(Theta)**2.)
    B = - np.log(Chi) + F_D(ReL,Theta) 
    return A/B



def F_L(ReL,Theta):#Theta en radians
    X = ReL*(1.-np.cos(Theta))
    Y = ReL*(1.+np.cos(Theta))
    CA1 = (2.-np.cos(Theta)+ np.cos(Theta)**2)*np.sin(Theta) / (2. * X)
    CA2 = (2.+np.cos(Theta)+ np.cos(Theta)**2)*np.sin(Theta) / (2. * Y)
    CA = CA1 * C(X) - CA2 * C(Y)
    CB = -3./2. +B(X)+B(Y)+1./2. * np.log(1-np.cos(Theta)**2) + gamma + np.log(ReL/4.)
    return CA * 1./np.sin(2*Theta) + CB

def L(Re,Chi,Theta):
    Re = np.array(Re)
    Chi = np.array(Chi)
    Theta = np.array(Theta)
    ReL = Re* Chi/2.
    try:
        if max(Theta) >= 3: 
            Theta = Theta / 360. *2. *math.pi
    except:
        if Theta >= 3:
            Theta = Theta / 360. *2. *math.pi
    A = 2. * math.pi *np.sin(Theta*2)
    B = -np.log(Chi) + F_L(ReL,Theta)
    return A/B

############ L et D normalse par la sphere de mm volume ##########

def L_norm(Re,Chi,Theta):
    return Chi**(2./3.)*(2./3.)**(1./3.)/(6.*math.pi) * L2(Re,Chi,Theta)

def D_norm(Re,Chi,Theta):
    return Chi**(2./3.)*(2./3.)**(1./3.)/(6.*math.pi) * D2(Re,Chi,Theta)

# Correcteur d'inertie

def Frond_D(ReL,Theta):   
    ReL = np.array(ReL)
    Theta = np.array(Theta)
    try:
        if max(Theta) >= 3: 
            Theta = Theta / 360. *2. *math.pi
    except:
        if Theta >= 3:
            Theta = Theta / 360. *2. *math.pi
    return F_D(ReL,Theta) - F_D(1e-7,Theta)

def Frond_L(ReL,Theta):   
    ReL = np.array(ReL)
    Theta = np.array(Theta)
    Theta = Theta / 360. *2. *math.pi
    return F_L(ReL,Theta) - F_L(1e-7,Theta)

# deuxieme forme pour D  et L

def D2(Re,Chi,Theta): 
    ReL = Chi*Re/2.
    return D(1e-7,Chi,Theta)*(1 + Frond_D(ReL,Theta) / np.log(Chi))

def L2(Re,Chi,Theta): # diverge en 0 
    ReL = Chi*Re/2.
    return L(1e-7,Chi,Theta)*(1 + Frond_L(ReL,Theta) / np.log(Chi))

## dans le repere perp te para 

def F_para_KC(Re,Chi,Theta):
    Re = np.array(Re)
    Chi = np.array(Chi)
    Theta = np.array(Theta)
    if max(Theta) >= 3: 
        Theta = Theta / 360. *2. *math.pi

    return np.cos(Theta) * D2(Re,Chi,Theta) +np.sin(Theta) *L2(Re,Chi,Theta) 

def F_perp_KC(Re,Chi,Theta):
    Re = np.array(Re)
    Chi = np.array(Chi)
    Theta = np.array(Theta)
    if max(Theta) >= 3: 
        Theta = Theta / 360. *2. *math.pi
    print 
    return +np.sin(Theta) * D2(Re,Chi,Theta) - np.cos(Theta) *L2(Re,Chi,Theta) 
