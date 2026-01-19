from sympy.tensor.tensor import TensorIndexType,tensor_indices, TensorHead, TensAdd,simplify_R, R,TensorSymmetry
from sympy.tensor.toperators import PartialDerivative
from sympy import *
from itertools import product
asym2 = TensorSymmetry.fully_symmetric(2)
asym3 = TensorSymmetry.fully_symmetric(3)

L = TensorIndexType("L",dim=3,dummy_name='L',metric_name='I')
I = L.delta
d = L.metric
x = TensorHead("x", [L])
idxs = tensor_indices("i j k l m k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 k_9 l_1 l_2 l_3 l_4 l_5",L)
i, j, k, l, m, k1, k2, k3, k4, k5, k6, k7, k8, k9, l1, l2, l3, l4, l5 = idxs


lamb, gam, E = symbols("L Gamma E")

Cf = symbols("C^f_{0:100}")
Cd = symbols("C^d_{0:100}")
const = [*Cf,*Cd]

def h(*idx):
    if idx:
        xs = [x(-i) for i in idx]
        return  PartialDerivative(R(-1),*xs)._perform_derivative(d,I)
    else:
        return R(-1)
    
def g(*idx):
    if idx:
        return simplify_R(R(len(idx)*2+1) * h(*idx))
    else:
        return 1
def sh(*idx):
    N = -Mul(*[-(3 +2*i) for i in range(len(idx))[:-1]])  
    if idx:
        return simplify_R(R(len(idx)+1) * h(*idx)) / N 
    else:
        return 1

    
def n(*idx):
    return x(*idx)*R(-1)

r = symbols("r")
Rto1 = [(R(n),1) for n in range(-20,20)]
Rto1 = [*Rto1,(r,1)]
Rtor = [(R(n),r**n) for n in range(-20,20)]
xtonr = [(x(idx),n(idx)*r) for idx in idxs]



def solve_all(Eqs):
    for eq in Eqs:
        display('1',eq)
    TEqs=[]
    for eq in Eqs:
        if isinstance(eq,TensAdd): 
            for e in eq.args:
                if e.args[0] == -1:
                    TEqs.append(e.args[0]*e.args[1])
                else:
                    TEqs.append(e.args[0])
        else:
            TEqs.append(eq.args[0])

    for e in TEqs:
        display(e)

    sols = solve(TEqs,const).items()
    for key,s  in sols:
        display(Eq(key, simplify(s)))

    return [(key, factor(s)) for key,s in sols]


def grad(f,*x):
    return PartialDerivative(f,*x)._perform_derivative(d,I)

def G(i,j):
    return (I(i,j)*R(-1)+x(i)*x(j)*R(-3))

#######

Pf2 = Cf[0]*h(k1)
Uf2 = (
    x(i) * Pf2 /2 
    + Cf[1]*I(i,k1)*h()
    + Cf[2]*h(i,k1)
)


Pd2 = Cd[0]*g(k1)
Ud2 = (
    x(i) * Pd2 /2 
    + Cd[1]*I(i,k1)*g()
    + Cd[2]*g(i,k1)
)     

########### SECOND ORDER 
Pf3 = Cf[1] * h(k1,k2)
Uf3 = (
    x(i) * Pf3 /2 
    + Cf[2]*h(k2)*I(k1,i)
    + Cf[3]*h(k1)*I(k2,i)
    + Cf[4]*h(i)*I(k1,k2)
    + Cf[5]*h(i,k1,k2)
)
Pd3 = Cd[1]*g(k1,k2)
Ud3 = (
    x(i) * Pd3 /2 
    + Cd[2]*g(k2)*I(k1,i)
    + Cd[3]*g(k1)*I(k2,i)
    + Cd[4]*g(i)*I(k1,k2)
    + Cd[5]*g(i,k1,k2)
)
## Consider a field u_i = U_ik1k2l K_k1k2l
# since d_k1l and d_k2l cancel the tensor doesnt contrain that
# It is also symmetric on k1k2 
Pf4 = (
    Cf[0]*I(k1,k2)*h(k3) 
    + Cf[1]*I(k3,k1)*h(k2) 
    + Cf[2]*I(k2,k3)*h(k1) 
    + Cf[3] * h(k1,k2,k3)
    )
Uf4 = (
    x(i) * Pf4 /2 
    + Cf[4]*I(k3,i)*I(k1,k2)*h()
    + Cf[5]*I(k1,i)*I(k3,k2)*h()
    + Cf[6]*I(k2,i)*I(k1,k3)*h()
    + Cf[7]*h(i,k1)*I(k2,k3)
    + Cf[8]*h(k3,i)*I(k1,k2)
    + Cf[9]*h(k2,i)*I(k3,k1)
    + Cf[10]*h(k1,k2)*I(k3,i)
    + Cf[11]*h(k1,k3)*I(i,k2)
    + Cf[12]*h(k2,k3)*I(i,k1) 
    + Cf[13]*h(i,k1,k2,k3)
)
Pd4 = (
    Cd[0]*I(k1,k2)*g(k3) 
    + Cd[1]*I(k3,k1)*g(k2) 
    + Cd[2]*I(k2,k3)*g(k1) 
    + Cd[3]*g(k1,k2,k3)
    )
Ud4 = (
    x(i) * Pd4 /2 
    + Cd[4]*I(k3,i)*I(k1,k2)*g()
    + Cd[5]*I(k1,i)*I(k3,k2)*g()
    + Cd[6]*I(k2,i)*I(k1,k3)*g()
    + Cd[7]*g(i,k1)*I(k2,k3)
    + Cd[8]*g(k3,i)*I(k1,k2)
    + Cd[9]*g(k2,i)*I(k3,k1)
    + Cd[10]*g(k1,k2)*I(k3,i)
    + Cd[11]*g(k1,k3)*I(i,k2)
    + Cd[12]*g(k2,k3)*I(i,k1) 
    + Cd[13]*g(i,k1,k2,k3)
)
## Consider a fiek3d u_i = P_k1k2k3m K_k1k2k3m and  u_i = U_ik1k2k3m K_k1k2k3m
# since d_k1m, d_k2m and d_k3m cancek3 the tensor doesnt contrain that
# It is ak3so symmetric on k1k2 
# for the pressure h2 dek3ta + h4  and this is it 
Pf5 = (
        +Cf[0]*I(k1,k2)*h(k3,k4)
        +Cf[1]*I(k1,k3)*h(k2,k4)
        +Cf[2]*I(k1,k4)*h(k3,k2)
        +Cf[3]*I(k3,k4)*h(k1,k2)
        +Cf[4]*I(k2,k4)*h(k1,k3)
        +Cf[5]*I(k3,k2)*h(k1,k4)
        +Cf[6]*h(k1,k2,k3,k4)
    )
## here we have I*I*h1 + I*H3 + h5
Uf5 = (
    x(i)*Pf5/2
    ### I*I have three perk4utation 
    + Cf[7]* I(i,k1)*I(k2,k3)*h(k4)
    + Cf[8]* I(i,k2)*I(k1,k3)*h(k4)
    + Cf[9]* I(i,k3)*I(k2,k1)*h(k4)
    ##
    + Cf[10]* I(i,k1)*I(k2,k4)*h(k3)
    + Cf[11]* I(i,k2)*I(k1,k4)*h(k3)
    + Cf[12]* I(i,k4)*I(k2,k1)*h(k3)
    ##
    + Cf[13]* I(i,k1)*I(k4,k3)*h(k2)
    + Cf[14]* I(i,k4)*I(k1,k3)*h(k2)
    + Cf[15]* I(i,k3)*I(k4,k1)*h(k2)
    ##
    + Cf[16]* I(i,k4)*I(k2,k3)*h(k1)
    + Cf[17]* I(i,k2)*I(k4,k3)*h(k1)
    + Cf[18]* I(i,k3)*I(k2,k4)*h(k1)
    ##
    + Cf[19]* I(k4,k1)*I(k2,k3)*h(i)
    + Cf[20]* I(k4,k2)*I(k1,k3)*h(i)
    + Cf[21]* I(k4,k3)*I(k2,k1)*h(i)
    ### I*h3 have 
    + Cf[22]* I(i,k1)*h(k2,k3,k4)
    + Cf[23]* I(k2,k1)*h(i,k3,k4)
    + Cf[24]* I(k3,k1)*h(k2,i,k4)
    + Cf[25]* I(k4,k1)*h(k2,k3,i)
    ##
    + Cf[26]* I(i,k2)*h(k1,k3,k4)
    + Cf[27]* I(k3,k2)*h(k1,i,k4)
    + Cf[28]* I(k4,k2)*h(k1,k3,i)
    ###
    + Cf[29]* I(i,k3)*h(k2,k1,k4)
    + Cf[30]* I(k4,k3)*h(k2,k1,i)
    ###
    + Cf[31]* I(i,k4)*h(k2,k3,k1)
    + Cf[32]* h(i,k1,k2,k3,k4)
) 
Pd5 = (
        +Cd[0]*I(k1,k2)*g(k3,k4)
        +Cd[1]*I(k1,k3)*g(k2,k4)
        +Cd[2]*I(k1,k4)*g(k3,k2)
        +Cd[3]*I(k3,k4)*g(k1,k2)
        +Cd[4]*I(k2,k4)*g(k1,k3)
        +Cd[5]*I(k3,k2)*g(k1,k4)
        +Cd[6]*g(k1,k2,k3,k4)
    )
## here we have I*I*h1 + I*H3 + h5
Ud5 = (
    x(i)*Pd5/2
    ### I*I have three perk4utation 
    + Cd[7]* I(i,k1)*I(k2,k3)*g(k4)
    + Cd[8]* I(i,k2)*I(k1,k3)*g(k4)
    + Cd[9]* I(i,k3)*I(k2,k1)*g(k4)
    ##
    + Cd[10]* I(i,k1)*I(k2,k4)*g(k3)
    + Cd[11]* I(i,k2)*I(k1,k4)*g(k3)
    + Cd[12]* I(i,k4)*I(k2,k1)*g(k3)
    ##
    + Cd[13]* I(i,k1)*I(k4,k3)*g(k2)
    + Cd[14]* I(i,k4)*I(k1,k3)*g(k2)
    + Cd[15]* I(i,k3)*I(k4,k1)*g(k2)
    ##
    + Cd[16]* I(i,k4)*I(k2,k3)*g(k1)
    + Cd[17]* I(i,k2)*I(k4,k3)*g(k1)
    + Cd[18]* I(i,k3)*I(k2,k4)*g(k1)
    ##
    + Cd[19]* I(k4,k1)*I(k2,k3)*g(i)
    + Cd[20]* I(k4,k2)*I(k1,k3)*g(i)
    + Cd[21]* I(k4,k3)*I(k2,k1)*g(i)
    ### I*h3 have 
    + Cd[22]* I(i,k1)*g(k2,k3,k4)
    + Cd[23]* I(k2,k1)*g(i,k3,k4)
    + Cd[24]* I(k3,k1)*g(k2,i,k4)
    + Cd[25]* I(k4,k1)*g(k2,k3,i)
    ##
    + Cd[26]* I(i,k2)*g(k1,k3,k4)
    + Cd[27]* I(k3,k2)*g(k1,i,k4)
    + Cd[28]* I(k4,k2)*g(k1,k3,i)
    ###
    + Cd[29]* I(i,k3)*g(k2,k1,k4)
    + Cd[30]* I(k4,k3)*g(k2,k1,i)
    ###
    + Cd[31]* I(i,k4)*g(k2,k3,k1)
    + Cd[32]* g(i,k1,k2,k3,k4)
) 





def compute_sol(X,sol):
    return simplify(X.canon_bp().subs(sol).doit())

def compute_surf(X):
    return X.contract_X((I,x,0)).subs(Rto1).doit().canon_bp()
        
def compute(X):
    return X.contract_X((I,x,0)).subs(Rtor).doit().canon_bp()

def compute_Sint(X):
    return simplify(compute_surf(X).int_X((I,x)).doit().contract_delta(I).canon_bp().doit())

sol4 = [(Cd[10], (lamb - 4)*(2*lamb + 3)/(12*(lamb + 1)*(lamb + 4))),
 (Cd[11], -(2*lamb**2 + 5*lamb + 8)/(12*(lamb + 1)*(lamb + 4))),
 (Cd[12], -(2*lamb**2 + 5*lamb + 8)/(12*(lamb + 1)*(lamb + 4))),
 (Cd[13], -1/(24*(lamb + 1))),
 (Cd[0], 0),
 (Cd[1], 0),
 (Cd[2], (2*lamb - 3)/(2*(lamb + 1))),
 (Cd[3], -1/(2*(lamb + 1))),
 (Cd[5], -1/(4*(lamb + 1))),
 (Cd[6], 0),
 (Cd[7], (lamb - 1)*(lamb + 5)/(6*(lamb + 1)*(lamb + 4))),
 (Cd[8], 0),
 (Cd[9], 0),
 (Cf[10], (13*lamb**2 + 10*lamb - 8)/(48*(lamb + 1)*(lamb + 4))),
 (Cf[11], -(lamb - 4)*(3*lamb + 2)/(48*(lamb + 1)*(lamb + 4))),
 (Cf[12], -(lamb - 4)*(3*lamb + 2)/(48*(lamb + 1)*(lamb + 4))),
 (Cf[13], lamb/(48*(lamb + 1))),
 (Cf[0], 0),
 (Cf[1], 0),
 (Cf[2], lamb/(4*(lamb + 1))),
 (Cf[3], (7*lamb + 2)/(24*(lamb + 1))),
 (Cf[4], 0),
 (Cf[5], -lamb/(8*(lamb + 1))),
 (Cf[6], 0),
 (Cf[7], (lamb - 1)/(6*(lamb + 1)*(lamb + 4))),
 (Cf[8], 0),
 (Cf[9], 0)]


sol2 = [(Cd[0], -5/(lamb + 1)),
 (Cd[1], -(2*lamb + 3)/(2*(lamb + 1))),
 (Cd[2], -1/(lamb + 1)),
 (Cf[0], (3*lamb + 2)/(2*(lamb + 1))),
 (Cf[1], -(3*lamb + 2)/(4*(lamb + 1))),
 (Cf[2], lamb/(4*(lamb + 1)))]

sol3 = [
(Cd[0],0),
(Cd[1], 7/(2*(lamb + 1))),
 (Cd[2], (2*lamb + 5)/(4*(lamb + 1))),
 (Cd[3], (2*lamb + 5)/(4*(lamb + 1))),
 (Cd[4], 0),
 (Cd[5], 5/(12*(lamb + 1))),
 (Cf[0], 0),
 (Cf[1], -(5*lamb + 2)/(3*(lamb + 1))),
 (Cf[2], 0),
 (Cf[3], 0),
 (Cf[4], 0),
 (Cf[5], -lamb/(6*(lamb + 1)))]

sol2B = [(Cd[0], 10/(3*(lamb + 1))),
 (Cd[1], 1/(3*(lamb + 1))),
 (Cd[2], 2/(3*(lamb + 1))),
 (Cf[0], 1/(3*(lamb + 1))),
 (Cf[1], -1/(6*(lamb + 1))),
 (Cf[2], 1/(6*(lamb + 1)))]

sol3B = [
(Cd[0], 0),
(Cd[1], -7/(10*(lamb + 1))),
 (Cd[2], (10*lamb + 1)/(60*(lamb + 1))),
 (Cd[3], -(10*lamb + 19)/(60*(lamb + 1))),
 (Cd[4], 1/(10*(lamb + 1))),
 (Cd[5], -1/(12*(lamb + 1))),
 (Cf[0],0),
 (Cf[1], -1/(5*(lamb + 1))),
 (Cf[2], 1/6),
 (Cf[3], -1/6),
 (Cf[4], 0),
 (Cf[5], -1/(30*(lamb + 1)))]

sol4B = [(Cd[0], -1/(6*(lamb + 1))),
 (Cd[10], (3*lamb + 5)/(42*(lamb + 1)*(lamb + 4))),
 (Cd[11], (3*lamb + 5)/(42*(lamb + 1)*(lamb + 4))),
 (Cd[12], -(2*lamb + 1)/(21*(lamb + 1)*(lamb + 4))),
 (Cd[13], 1/(210*(lamb + 1))),
 (Cd[1], -1/(6*(lamb + 1))),
 (Cd[2], 2/(3*(lamb + 1))),
 (Cd[3], 2/(35*(lamb + 1))),
 (Cd[4], -1/(60*(lamb + 1))),
 (Cd[5], 1/(15*(lamb + 1))),
 (Cd[6], -1/(60*(lamb + 1))),
 (Cd[7], (3*lamb + 19)/(42*(lamb + 1)*(lamb + 4))),
 (Cd[8], -(lamb + 11)/(84*(lamb + 1)*(lamb + 4))),
 (Cd[9], -(lamb + 11)/(84*(lamb + 1)*(lamb + 4))),
 (Cf[0], -1/(60*(lamb + 1))),
 (Cf[10], (5*lamb + 6)/(84*(lamb + 1)*(lamb + 4))),
 (Cf[11], (5*lamb + 6)/(84*(lamb + 1)*(lamb + 4))),
 (Cf[12], -(9*lamb + 8)/(84*(lamb + 1)*(lamb + 4))),
 (Cf[13], 1/(420*(lamb + 1))),
 (Cf[1], -1/(60*(lamb + 1))),
 (Cf[2], 1/(15*(lamb + 1))),
 (Cf[3], 1/(42*(lamb + 1))),
 (Cf[4], 1/(120*(lamb + 1))),
 (Cf[5], -1/(30*(lamb + 1))),
 (Cf[6], 1/(120*(lamb + 1))),
 (Cf[7], -(lamb - 3)/(42*(lamb + 1)*(lamb + 4))),
 (Cf[8], (3*lamb - 2)/(168*(lamb + 1)*(lamb + 4))),
 (Cf[9], (3*lamb - 2)/(168*(lamb + 1)*(lamb + 4)))]

sol4Bsym =[(Cd[0], 1/(9*(lamb + 1))),
 (Cd[10], 1/(63*(lamb + 1))),
 (Cd[11], 1/(63*(lamb + 1))),
 (Cd[12], 1/(63*(lamb + 1))),
 (Cd[13], 1/(210*(lamb + 1))),
 (Cd[1], 1/(9*(lamb + 1))),
 (Cd[2], 1/(9*(lamb + 1))),
 (Cd[3], 2/(35*(lamb + 1))),
 (Cd[4], 1/(90*(lamb + 1))),
 (Cd[5], 1/(90*(lamb + 1))),
 (Cd[6], 1/(90*(lamb + 1))),
 (Cd[7], 1/(63*(lamb + 1))),
 (Cd[8], 1/(63*(lamb + 1))),
 (Cd[9], 1/(63*(lamb + 1))),
 (Cf[0], 1/(90*(lamb + 1))),
 (Cf[10], 1/(252*(lamb + 1))),
 (Cf[11], 1/(252*(lamb + 1))),
 (Cf[12], 1/(252*(lamb + 1))),
 (Cf[13], 1/(420*(lamb + 1))),
 (Cf[1], 1/(90*(lamb + 1))),
 (Cf[2], 1/(90*(lamb + 1))),
 (Cf[3], 1/(42*(lamb + 1))),
 (Cf[4], -1/(180*(lamb + 1))),
 (Cf[5], -1/(180*(lamb + 1))),
 (Cf[6], -1/(180*(lamb + 1))),
 (Cf[7], 1/(252*(lamb + 1))),
 (Cf[8], 1/(252*(lamb + 1))),
 (Cf[9], 1/(252*(lamb + 1)))]