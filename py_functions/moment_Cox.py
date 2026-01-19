import numpy as np
import scipy.special as sc
import math
# Couple proposee par Khayat  Cox
gamma = 0.5772156649

def T_simple(Re,Chi,Theta):
    Theta = np.array(Theta) /360 * 2*math.pi
    T = 5.*math.pi/12. *Chi*Re / np.log(Chi)**2 * np.sin(2*Theta)
    return T

def T_ad(Re,Chi,Theta):
    Theta = np.array(Theta) /360 * 2*math.pi
    Re = np.array(Re)
    Chi = np.array(Chi)
    A = 5.*math.pi/(12.*4.) * Re/(1.+Chi*Re**1.1)**0.5
    B = Chi/np.log(Chi)**2
    C = (13.5-30.*Re**(1./2.))/(Chi*np.log(Chi)**3)*np.exp(-0.7*Re)
    return A * (B+C) * np.sin(2*Theta)
def T_adchi2(Re,Chi,Theta):
    Theta = np.array(Theta) /360 * 2*math.pi
    Re = np.array(Re)
    Chi = np.array(Chi)
    A = 5.*math.pi/(12.*4.) * Re/(.5 +Chi*Re**1.1)**0.5
    B = Chi/np.log(2*Chi)**2
    C = (Chi*2*np.log(2)+(13.5-30.*Re**(1./2.))*np.exp(-0.7*Re))/(Chi*np.log(2*Chi)**3)
    return A * (B+C) * np.sin(2*Theta)

def Zg(x):
    x = np.array(x)
    Z = 2./x * (1.+(np.exp(-x)-1.)/x)
    return Z
def Cg(X):
    return (sc.expn(1,X)+np.log(X)+gamma)/X

def F_g(ReL,Theta):
    Theta = np.array(Theta)
    Theta = Theta /360. * 2.*math.pi
    X =ReL*(1.-np.cos(Theta))
    Y =ReL*(1.+np.cos(Theta))
    B = Zg(X) -Cg(X) + Zg(Y) -Cg(Y)
    return (np.cos(Theta) * B + Zg(Y) - Zg(X)) * np.sin(Theta)

# Torque KC
def T_CC(Re,Chi,Theta) :
    Re = np.array(Re)
    Chi = np.array(Chi)
    ReL = Re * Chi / 2.
    A = - 0.5* math.pi/np.log(Chi)**2.
    return  A * F_g(ReL,Theta)


def T_CC_lim(Re,Chi,Theta) :
    Theta = np.array(Theta) /360. * 2.*math.pi
    Re = np.array(Re)
    Chi = np.array(Chi)
    ReL = Re * Chi / 2.
    A = 5./24.*math.pi * np.sin(2*Theta) *ReL
    return  A /np.log(Chi)**2.

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
def F_CC_perunitoflength(Re):
    gamma = 0.5772156649
    A = - 8. * math.pi * Re**(-1) 
    B = np.log(1./2. * Re)+ gamma -0.5 -2.*np.log(2.)
    return A/B


def C_d_ref(Re):
    Re = np.array(Re)
    return 9.689*Re**(-0.78)*(1+0.147*Re**(0.82))


def F_90_myfit(Re,Chi):
    C3,C4,C6,C7,C9,C10,C11,C12 =2.83304331e+02,  1.97993497e+01, -1.48065106e+02,  4.37651735e+00,3.00000000e+00,  1.00000000e+00,  4.21147595e-01,  1.82966549e-01 
    Re = np.array(Re)
    Chi = np.array(Chi)
    Re_L = Chi*Re/2.
    f_3 = C3 *np.exp(-C4/Re_L**C11)
    f_4 = C6 *np.exp(-C7/Re_L**C12)
    # f_3 = C3 *np.exp(-C4/Chi)*Re**C12
    # f_4 = C6 *np.exp(-C7/Chi)*Re**C11
    # f_4 = C6*(Chi* Re)**(C3*(Chi)**C4)
    # f_4 = C6*(Chi* Re)**3.
    # f_3 = C3*(Chi* Re)**2.

    f_perp =sc.expn(1,Re_L)+np.log(Re_L)-(np.exp(-Re_L)-1.)/(Re_L)+gamma-1.
    A = 1./np.log(2*Chi)
    B =( - 0.19315 +f_perp) /np.log(2*Chi)**2
    C = (0.21484 + f_3 +2*f_perp*np.log(2) )/np.log(2*Chi)**3
    # C = (0.21484 + f_3 )/np.log(2*Chi)**3
    D = (0.38735 + 3*np.log(2.)*f_3 +3*f_perp*np.log(2)**2. +f_4)/np.log(2*Chi)**4
    # D = (0.38735  +f_3+f_4)/np.log(2*Chi)**4
    # return 4./3*(Chi)/Chi**(1./3)*(2./3)**(1./3)*(A+B+chicoef*(C+D))
    return 4.*math.pi*((A+B+C+D)-0.56819516/(Chi-0.5)**1.75*Chi**(-2./3.))
def F_90_KC2(Re,Chi):
    C3,C4,C6,C7,C9,C10,C11,C12 =2.83304331e+02,  1.97993497e+01, -1.48065106e+02,  4.37651735e+00,3.00000000e+00,  1.00000000e+00,  4.21147595e-01,  1.82966549e-01 
    Re = np.array(Re)
    Chi = np.array(Chi)
    Re_L = Chi*Re/2.
    f_3 = C3 *np.exp(-C4/Re_L**C11)
    f_4 = C6 *np.exp(-C7/Re_L**C12)
    # f_3 = C3 *np.exp(-C4/Chi)*Re**C12
    # f_4 = C6 *np.exp(-C7/Chi)*Re**C11
    # f_4 = C6*(Chi* Re)**(C3*(Chi)**C4)
    # f_4 = C6*(Chi* Re)**3.
    # f_3 = C3*(Chi* Re)**2.

    f_perp =sc.expn(1,Re_L)+np.log(Re_L)-(np.exp(-Re_L)-1.)/(Re_L)+gamma-1.
    A = 1./np.log(2*Chi)
    B =( - 0.19315 +f_perp) /np.log(2*Chi)**2
    C = (0.21484 + f_3 +2*f_perp*np.log(2) )/np.log(2*Chi)**3
    # C = (0.21484 + f_3 )/np.log(2*Chi)**3
    D = (0.38735)/np.log(2*Chi)**4
    # D = (0.38735  +f_3+f_4)/np.log(2*Chi)**4
    # return 4./3*(Chi)/Chi**(1./3)*(2./3)**(1./3)*(A+B+chicoef*(C+D))
    return 4.*math.pi*((A+B+C+D)-0.56819516/(Chi-0.5)**1.75*Chi**(-2./3.))


def F_90_chi(Chi):
    Chi = np.array(Chi)
    A = 1./np.log(2*Chi)
    B = - 0.19315 /np.log(2*Chi)**2
    C = 0.21484/np.log(2*Chi)**3
    D = 0.38735 /np.log(2*Chi)**4
    return 4./3*(Chi)/Chi**(1./3)*(2./3)**(1./3)*(A+B+C+D)


def F_90_chifit(Chi):
    Chi = np.array(Chi)
    A = 1./np.log(2*Chi)
    B = - 0.19315 /np.log(2*Chi)**2
    C = 0.21484/np.log(2*Chi)**3
    D = 0.38735 /np.log(2*Chi)**4
    return 4*math.pi*((A+B+C+D)-0.56819516/(Chi-0.5)**1.75*Chi**(-2./3.))

def F_90_KC(Re,Chi):
    # def F_90_myfit(X,C1,C2,C3,C4,C5,C6):
    Re = np.array(Re)
    Chi = np.array(Chi)
    Re_L =Re*Chi/2.
    f_perp =sc.expn(1,Re_L)+np.log(Re_L)-(np.exp(-Re_L)-1.)/(Re_L)+gamma-1.
    A = 1./np.log(2*Chi)
    B =( - 0.19315 +f_perp) /np.log(2*Chi)**2
    C = (0.21484)/np.log(2*Chi)**3
    D = (0.38735)/np.log(2*Chi)**4
    # return 4./3*(Chi)/Chi**(1./3)*(2./3)**(1./3)*(A+B+chicoef*(C+D))
    return 4.*math.pi*((A+B+C+D)-0.56819516/(Chi-0.5)**1.75*Chi**(-2./3.))

def F_0_KC(Re,Chi):
    Re = np.array(Re)
    Chi = np.array(Chi)
    f_A = (sc.expn(1,Chi*Re) + np.log(Chi*Re) - np.exp(-Chi*Re) + 1 +gamma)/(Chi*Re)
    f_B = sc.expn(1,Chi*Re) + np.log(Chi*Re) + gamma -2.
    f_para = 1./2. * (f_A + f_B)
    A = 1./np.log(Chi)
    B = (0.80685 + f_para -np.log(2.))* 1./np.log(Chi)**2
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B)#-1.4/(Chi-0.5)**1.75
    return F_para2


def F2A(Re,Chi):
    ReL=Re*Chi/2
    # return 0.0869545 * ReL
    return 0.110207 * ReL # normalisation

def F2B(Re,Chi):
    ReL=Re*Chi/2
    return -10**(-1.41392 + 0.508603 *np.log(ReL) - 0.0396166* np.log(ReL)**2)
def F2B(Re,Chi):
    ReL=Re*Chi/2
    C1,C2,C3=1.2603546765039138, 0.815606561137388, 2.34131761572026
    F2B = 0.215674107599209*np.log(1+(ReL)**(C1))/((ReL)**(C1-2)+C3*(ReL)**(C2-1))
    return -F2B#+C4*(Re)**(C2))

def F2C():
    return 0.82854 + np.log(2)**2-0.80685*np.log(2)*2

# exp en ln Chi
def F_0_KC3(Re,Chi):
    Re = np.array(Re)
    Chi = np.array(Chi)
    f_A = (sc.expn(1,Chi*Re) + np.log(Chi*Re) - np.exp(-Chi*Re) + 1 +gamma)/(Chi*Re)
    f_B = sc.expn(1,Chi*Re) + np.log(Chi*Re) + gamma -2.
    f_para = 1./2. * (f_A + f_B)
    A = 1./np.log(Chi)
    B = (0.80685 + f_para -np.log(2.))* 1./np.log(Chi)**2
    C = (F2A(Re,Chi)+F2B(Re,Chi)+F2C())* 1./np.log(Chi)**3
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C)#-1.4/(Chi-0.5)**1.75
    return F_para2

# exp en ln 2*Chi
def F_0_KC2_2C(Re,Chi):
    Re = np.array(Re)
    Chi = np.array(Chi)
    f_A = (sc.expn(1,Chi*Re) + np.log(Chi*Re) - np.exp(-Chi*Re) + 1 +gamma)/(Chi*Re)
    f_B = sc.expn(1,Chi*Re) + np.log(Chi*Re) + gamma -2.
    f_para = 1./2. * (f_A + f_B)
    A = 1./np.log(2*Chi)
    B = (0.80685 + f_para)* 1./np.log(2*Chi)**2
    C = 0.82854 /np.log(2*Chi)**3
    D = (1.45243) *1./np.log(2*Chi)**4
    F_para2 = 2.*math.pi*(A+B+C+D-2.34270029/(Chi-0.5)**1.75*Chi**(-2./3.))
    return F_para2
# exp en ln 2*Chi
def F_0_KC3_2C(Re,Chi):
    Re = np.array(Re)
    Chi = np.array(Chi)
    f_A = (sc.expn(1,Chi*Re) + np.log(Chi*Re) - np.exp(-Chi*Re) + 1 +gamma)/(Chi*Re)
    f_B = sc.expn(1,Chi*Re) + np.log(Chi*Re) + gamma -2.
    f_para = 1./2. * (f_A + f_B)
    A = 1./np.log(2*Chi)
    B = (0.80685 + f_para)* 1./np.log(2*Chi)**2
    C = (F2A(Re,Chi)+F2B(Re,Chi) +0.82854 +f_para*np.log(4) )* 1./np.log(2*Chi)**3
    D = (1.45243) *1./np.log(2*Chi)**4
    F_para2 = 2.*math.pi*(A+B+C+D-2.34270029/(Chi-0.5)**1.75*Chi**(-2./3.))
    return F_para2
    # (3*f_para*np.log(2)**2   +3*np.log(2)* (F2A(Re,Chi)+F2B(Re,Chi)))*

    
def F_0(Re,Chi):
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





# sans inertie avec fit en Chi
def F_0_modi(Chi):
    Chi = np.array(Chi)
    f_3 = 1#(Chi*Re)**(0.07*Chi**0.5)
    f_4 = 1#(Chi* Re)**(0.03*Chi**0.9)
    A = 1./np.log(2*Chi)
    B = 0.80685* 1./np.log(2*Chi)**2
    C = 0.82854* 1./np.log(2*Chi)**3
    D = 1.45243* 1./np.log(2*Chi)**4
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D)-1.4/(Chi-0.5)**1.75
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D)-1.4/(Chi-0.5)**1.75
    return F_para2

# 4th order sans correction 
def F_04th(Chi):
    Chi = np.array(Chi)
    A = 1./np.log(2*Chi)
    B = 0.80685* 1./np.log(2*Chi)**2
    C = 0.82854* 1./np.log(2*Chi)**3
    D = 1.45243* 1./np.log(2*Chi)**4
    F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D)
    return F_para2

def F_04th2(Chi):
    a1,a2,a3,a4 = 1.,0.80685,0.82854,1.45243
    A = a1/np.log(Chi)
    B = (a2-np.log(2.)*a1)/np.log(Chi)**2
    C = (a3+np.log(2.)**2*a1-2*np.log(2.)*a2)/np.log(Chi)**3
    D = (a4-np.log(2.)**3*a1+  3*np.log(2.)**2*a2 - 3*np.log(2.)*a3)/np.log(Chi)**4
    return 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D) 

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
    return F_90_myfit(Re,Chi) * np.sin(theta)

def F_perp_Modi(Re,Chi,theta):
    theta = np.array(theta) /360 * 2*math.pi
    return (F_90_myfit(Re,Chi) + 1./16. *np.cos(theta)**2  *Chi*Re / (np.log(Chi)**2)) * np.sin(theta)

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

def F_Ds(Theta):#theta en radians
    A = 2*(1-4*np.log(2)) + (1+4*np.log(2.))*np.cos(Theta)**2
    return A/(2*(2-np.cos(Theta)))

def F_Ls(Theta):#theta en radians
    return -1./2.-2*np.log(2.)


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

    A = - 2. * math.pi * (2. - np.cos(Theta)**2.)
    B = - np.log(Chi) + F_D(ReL,Theta) 
    return A/B

def Ds(Chi,Theta): # diverge en 0
    Theta = np.array(Theta)
    Chi = np.array(Chi)
    try:
        if max(Theta) >= 3: 
            Theta = Theta / 360. *2. *math.pi
    except:
        if Theta >= 3:
            Theta = Theta / 360. *2. *math.pi

    A = - 2. * math.pi * (2. - np.cos(Theta)**2.)
    B = - np.log(Chi) + F_Ds(Theta) 
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
    A = 1. * math.pi *np.sin(Theta*2)
    B = -np.log(Chi) + F_L(ReL,Theta)
    return A/B

def Ls(Chi,Theta):
    Chi = np.array(Chi)
    Theta = np.array(Theta)
    try:
        if max(Theta) >= 3: 
            Theta = Theta / 360. *2. *math.pi
    except:
        if Theta >= 3:
            Theta = Theta / 360. *2. *math.pi
    A = 1. * math.pi *np.sin(Theta*2)
    B = -np.log(Chi) + F_Ls(Theta)
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
    return F_D(ReL,Theta) - F_Ds(Theta)

def Frond_L(ReL,Theta):   
    ReL = np.array(ReL)
    Theta = np.array(Theta)
    try:
        if max(Theta) >= 3: 
            Theta = Theta / 360. *2. *math.pi
    except:
        if Theta >= 3:
            Theta = Theta / 360. *2. *math.pi
    return F_L(ReL,Theta) - F_Ls(Theta)

# deuxieme forme pour D  et L
# Drag and lift of K and C 
def D2(Re,Chi,Theta): 
    ReL = Chi*Re/2.
    return Ds(Chi,Theta)*(1 + Frond_D(ReL,Theta) / np.log(Chi))

def L2(Re,Chi,Theta): # diverge en 0 
    ReL = Chi*Re/2.
    return Ls(Chi,Theta) *(1 + Frond_L(ReL,Theta) / np.log(Chi))



# 
def D2_old(Re,Chi,Theta): 
    ReL = Chi*Re/2.
    return D(1e-7,Chi,Theta)*(1 + Frond_D(ReL,Theta) / np.log(Chi))

def L2_old(Re,Chi,Theta): # diverge en 0 
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
    return - np.sin(Theta) * D2(Re,Chi,Theta) + np.cos(Theta) * L2(Re,Chi,Theta) 

def T_myfit(Re,Chi):#
    C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11 = 0.57679654,  0.48551775, -0.05001353,  0.25083972,  0.41980062,  0.39046998,1.27116393,  0.53993132,  0.59604704,  1. ,         1.62802176
    Re = np.array(Re)
    Chi = np.array(Chi)

    
    # R = C1*(Re)**C2+C3
    R = C1*(Re)**(C2)+C3
    R2 = C4*np.sin((Chi + 30./2. * np.tanh( (np.log(Re /C5)))  ) * math.pi /30.) 
    R3 = -C6* np.exp(-C7*( np.log(Re/C8))**2)
    R4 =Re**C9/np.log(Chi*C11)
    

    return  R * (1+R2*(1+R3) *(1+R4))


def f_0(Re,Chi):
    return 0.01211161*Chi*Re

def f_90(Re,Chi):
    return 0.01251501*Chi*Re

import scipy.special as sp
def T_perp(Reo,Chi):
    Reo = np.array(Reo)
    Chi = np.array(Chi)
    eps = 1 /np.log(2*Chi)
    f = 0.0337*Chi**(1.3)*Reo**(0.9)
    A = (11./6. -math.log(2) + f)
    B = (161./36. -math.pi**2/12.-11./3.*np.log(2.)+np.log(2.)**2)
    E = -5./4. * sp.zeta(3.) + 1033./72.-np.log(2.)**3 +11./2. *np.log(2.)**2 - 161./12. * np.log(2.) - math.pi**2 *(11./24. -1./4.*np.log(2.))
    D = (1-1/(2.*Chi)**(1.2))**5*E
    T = eps +eps**2*A + eps**3 *B +eps**4*D
    return T

def F_0_myfit(Re,Chi):
    C1,C2,C3,C4,C5,C6 = 0.06      ,  0.5       ,  0.76241099,  0.9       ,  1.        ,-0.63645881
    import numpy as np
    import math
    Re = np.array(Re)
    Chi = np.array(Chi)
    f_A = (sc.expn(1,Chi*Re) + np.log(Chi*Re) - np.exp(-Chi*Re) + 1 +gamma)/(Chi*Re)
    f_B = sc.expn(1,Chi*Re) + np.log(Chi*Re) + gamma -2.
    f_para = 1./2. * (f_A + f_B)
    # (3*f_para*np.log(2)**2   +3*np.log(2)* (F2A(Re,Chi)+F2B(Re,Chi)))*
    f_3 = C5*(Chi*Re/2.)**(C1*Chi**C2)
    # f_4 = C6*(Chi* Re/2.)**(C3*(Chi)**C4) * Chi**(-1.)
    f_4 = C6*(Chi* Re/2.)**(C3)
    # f_4 = C6*(Chi* Re)**(3)*Chi**C3
    A = 1./np.log(2*Chi)
    B = (0.80685 + f_para )* 1./np.log(2*Chi)**2
    C = (0.82854 +( F2A(Re,Chi)+F2B(Re,Chi) +2*f_para*np.log(2)))* 1./np.log(2*Chi)**3
    D = (1.45243  + (3*np.log(2.)*( F2A(Re,Chi)+F2B(Re,Chi)) +3*f_para*np.log(2)**2.)+f_4)* 1./np.log(2*Chi)**4
    # F_para2 = 1./3*(Chi)**(2./3)*(16./3)**(1./3)*(A+B+C+D)-1.4/(Chi-0.5)**1.75
    F_para2 = 2.*math.pi*((A+B+C+D)-2.34270029/(Chi-0.5)**1.75*Chi**(-2./3.))
    return F_para2
def F_90_myfit(Re,Chi):
    C3,C4,C6,C7,C9,C10,C11,C12 = 2.83304331e+02,  1.97993497e+01, -1.48065106e+02,  4.37651735e+00,3.00000000e+00,  1.00000000e+00,  4.21147595e-01,  1.82966549e-01
    Re = np.array(Re)
    Chi = np.array(Chi)
    Re_L = Chi*Re/2.
    f_3 = C3 *np.exp(-C4/Re_L**C11)
    f_4 = C6 *np.exp(-C7/Re_L**C12)
    # f_3 = C3 *np.exp(-C4/Chi)*Re**C12
    # f_4 = C6 *np.exp(-C7/Chi)*Re**C11
    # f_4 = C6*(Chi* Re)**(C3*(Chi)**C4)
    # f_4 = C6*(Chi* Re)**3.
    # f_3 = C3*(Chi* Re)**2.

    f_perp =sc.expn(1,Re_L)+np.log(Re_L)-(np.exp(-Re_L)-1.)/(Re_L)+gamma-1.
    A = 1./np.log(2*Chi)
    B =( - 0.19315 +f_perp) /np.log(2*Chi)**2
    C = (0.21484 + f_3 +2*f_perp*np.log(2) )/np.log(2*Chi)**3
    # C = (0.21484 + f_3 )/np.log(2*Chi)**3
    D = (0.38735 + 3*np.log(2.)*f_3 +3*f_perp*np.log(2)**2. +f_4)/np.log(2*Chi)**4
    # D = (0.38735  +f_3+f_4)/np.log(2*Chi)**4
    # return 4./3*(Chi)/Chi**(1./3)*(2./3)**(1./3)*(A+B+chicoef*(C+D))
    return 4.*math.pi*((A+B+C+D)-0.56819516/(Chi-0.5)**1.75*Chi**(-2./3.))
def F_perp_modi_myfit(Re,Chi,theta):
    C1,C2 = 0.0215796 , 0.68514292
    theta = np.array(theta) /360 * 2*math.pi
    return (F_90_myfit(Re,Chi) -2*math.pi* C1 *np.cos(theta)**2  *(Chi*Re/2.)**C2) * np.sin(theta)
def F_para_modi_myfit(Re,Chi,theta):
    C1,C2 = 0.03978024, 0.54469737
    theta = np.array(theta) /360 * 2*math.pi
    # return (F_0_myfit(Re,Chi) + C1 *np.sin(theta)**2  *Re*Chi )* np.cos(theta)
    return (F_0_myfit(Re,Chi) +2*math.pi* C1 *np.sin(theta)**2  *(Re*Chi/2.)**C2)* np.cos(theta)


def T_myfit(Re,Chi):
    Re = np.array(Re)
    Chi = np.array(Chi)
    ReL = Re*Chi/2.
    CA = 3.
    A = 5./24.*math.pi *ReL/(1.+ ReL**1.99138079 )**0.33120667
    B = 1.
    C = 2.24375843 -1.81333523*ReL**0.54383989 
    D = -3.60337783 + 8.85404603*ReL**0.53767024
    E = -14.30056829*(ReL/Chi)**0.4484439
    return A * (B/np.log(CA*Chi)**2+C /np.log(CA*Chi)**3+D /np.log(CA*Chi)**4 + E/np.log(CA*Chi)**5)

import numpy.linalg as LA
def R_FU(Re: float,Chi: float, p: np.array,u: np.array):
    un = u/LA.norm(u,2)
    pp = np.outer(p,p)
    ipp = np.identity(3)  -  np.outer(p,p)
    R_linear = F_0_myfit(Re,Chi) * pp + F_90_myfit(Re,Chi) * ipp 
    R_inertie_para = f_0(Re,Chi) * np.dot(un,np.dot(ipp,un)) * pp
    R_inertie_perp = f_90(Re,Chi)* np.dot(un,np.dot(pp,un)) * ipp
    return(R_linear + R_inertie_para + R_inertie_perp)

def R_TU(Re: float,Chi: float, p: np.array,u: np.array):
    un = u/LA.norm(u,2)
    return  2* np.outer(np.cross(p,un) , p)  * T_myfit(Re,Chi)

def R_Tomega(Reo: float,Chi: float,p: np.array):
    ipp = np.identity(3)  -  np.outer(p,p)
    return  (T_perp(Reo,Chi) * ipp )



 








# line force 

def f_d(s:float,Chi: float,Re:float,Theta:float):
    Re = np.array(Re)
    Chi = np.array(Chi)
    Theta = np.array(Theta) / 180. *math.pi
    Re_L = Chi*Re/2.
    Beta = np.cos(Theta)
    e = 1
    ct = np.cos(Theta)
    A = ct*Beta - 2* e
    X = 1./2. * Re_L * (1-ct)*(1+s)
    Y = 1./2. * Re_L * (1+ct)*(1-s)
    B1 =  1./4.*(2*ct*e - (2-ct+ct**2)*Beta)*((1.- np.exp(-X))/X - 1)
    B2 = -1./4.*(2*ct*e - (2+ct+ct**2)*Beta)*((1.-np.exp(-Y))/Y - 1)
    B3 = -1./2. * (ct * Beta - 2*e) * (sc.expn(1,X) + np.log(1-ct))
    B4 = -1./2. * (ct * Beta - 2*e) * (sc.expn(1,Y) + np.log(1+ct))
    B5 = -(ct* Beta - 2 * e) * (gamma + np.log(Re_L/4.)) + 3./2. * ct * Beta - e
    B = B1+B2+B3+B4+B5
    return math.pi*2*(-A/np.log(Chi) + B/np.log(Chi)**2)

def f_l(s:float,Chi: float,Re:float,Theta:float):
    Re = np.array(Re)
    Chi = np.array(Chi)
    Theta = np.array(Theta) / 180. *math.pi
    Re_L = Chi*Re/2.
    Beta = np.sin(Theta)
    e = 0.
    ct = np.cos(Theta)
    A = ct*Beta - 2* e
    X = 1./2. * Re_L * (1-ct)*(1+s)
    Y = 1./2. * Re_L * (1+ct)*(1-s)
    B1 =  1./4.*(2*ct*e - (2-ct+ct**2)*Beta)*((1.- np.exp(-X))/X - 1)
    B2 = -1./4.*(2*ct*e - (2+ct+ct**2)*Beta)*((1.-np.exp(-Y))/Y - 1)
    B3 = -1./2. * (ct * Beta - 2*e) * (sc.expn(1,X) + np.log(1-ct))
    B4 = -1./2. * (ct * Beta - 2*e) * (sc.expn(1,Y) + np.log(1+ct))
    B5 = -(ct* Beta - 2 * e) * (gamma + np.log(Re_L/4.)) + 3./2. * ct * Beta - e
    B = B1+B2+B3+B4+B5
    return math.pi*2*(-A/np.log(Chi) + B/np.log(Chi)**2)


def m(s:float,Chi: float,Re:float,Theta:float):
    ct = np.cos(Theta/180*math.pi)
    st = np.sin(Theta/180*math.pi)
    return s*(ct*f_l(s,Chi,Re,Theta) - st *f_d(s,Chi,Re,Theta))/2.

def T_chawng(Chi):
    return -(1. - 1./2.*Chi**(-1.) + 1./2.*(Chi)**(-2.)) / 2.
    
    
    

# def FparaKC(Re: float,Chi: float,Theta: float,F0: function,para: tuple):
#     Theta = np.array(Theta) /360 * 2*math.pi
#     ReL = Chi*Re/2.
#     A = - 58**np.cos(Theta) + 25.*np.cos(3*Theta)+ np.cos(5*Theta)+ 202*np.sin(Theta)- 23*np.sin(3* Theta)- np.sin(5*Theta)
#     B = 64.*(np.cos(2*Theta) - 3)
#     return 1# (F0((Re,Chi),*para) + A/B * ReL /np.log(Chi)**2 )