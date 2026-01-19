from sympy.tensor.tensor import TensorIndexType,tensor_indices, TensorHead, TensAdd, Tensor,TensMul,TensExpr, R,TensorSymmetry, contract_X
from sympy.tensor.toperators import PartialDerivative
from sympy import pi
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
phi = Symbol('phi', positive=True)

lamb = Symbol('L', positive=True)
Cf = symbols("C^f_{0:100}")
Cd = symbols("C^d_{0:100}")
const = [*Cf,*Cd]

# def simplify_R(expr):
#     if expr.args:
#         if isinstance(expr,TensAdd):
#             return TensAdd(*[simplify_R(x) for x in expr.args])
#         if isinstance(expr,Add):
#             return Add(*[simplify_R(x) for x in expr.args])
            
#     if len([arg for arg in expr.args if  isinstance(arg,R)]) > 1:
#         exprSr = Mul(*[arg for arg in expr.args if not isinstance(arg,R)])
#         exprWr = [arg for arg in expr.args if  isinstance(arg,R)][0]
#         for exp in [arg for arg in expr.args if  isinstance(arg,R)][1:]:
#             exprWr = exprWr*exp
            
#         return exprSr*exprWr
        
#     if isinstance(expr,TensMul):
#         return TensMul(*[simplify_R(exp) if isinstance(exp,(Add,TensAdd,Pow)) else exp for exp in expr.args])
    
#     if isinstance(expr,Mul):
#         return Mul(*[simplify_R(exp) if isinstance(exp,(Add,Pow)) else exp for exp in expr.args])
    
#     if isinstance(expr,Pow) and isinstance(expr.args[0],R):
#         if type(expr.args[0])==R:
#             rN = expr.args[0].N
#             rXo = expr.args[0].Xo
#             n = expr.args[1]
#             return R(rN*n,rXo)
#         else:
#             return expr
#     return expr.doit()

def grad(f,*x):
    return PartialDerivative(f,*x)._perform_derivative(d,I)
## solution with zero pressure 
def UH(*idx):
    i = idx[0]
    idx = idx[1:]
    if idx:
        xs = [x(-i) for i in idx]
        return  - grad(R(-1),x(-i),*xs) /4
    else:
        return  - grad(R(-1),x(-i)) /4
def SIGH(*idx):
    i = idx[0]
    j = idx[1]
    idx = idx[2:]
    if idx:
        xs = [x(-i) for i in idx]
        return  - 2*grad(R(-1),x(-i),x(-j),*xs) /4
    else:
        return  - 2*grad(R(-1),x(-i),x(-j)) /4
#greenn func of order N
def UF(*idx):
    i = idx[0]
    j = idx[1]
    idx = idx[2:]
    if idx:
        xs = [x(-i) for i in idx]
        return  (grad(R(1),x(-i),x(-j),*xs) -grad(R(1),x(-k1),x(k1),*xs)*I(i,j))/8
    else:
        return   (grad(R(1),x(-i),x(-j)) -grad(R(1),x(-k1),x(k1))*I(i,j))/8
##Singular stresslet 
def SIGF(*idx):
    i = idx[0]
    j = idx[1]
    k1 = idx[2]
    idx = idx[3:]
    if idx:
        xs = [x(-i) for i in idx]
        return  (2*grad(R(1) ,x(-k1),x(-i),x(-j), *xs) 
                 -I(k1,i)*grad(R(1) ,x(-l1),x(l1),x(-j) ,*xs) 
                 -I(i,j)*grad(R(1)  ,x(-k1),x(-l1),x(l1),*xs) 
                 -I(k1,j)*grad(R(1) ,x(-l1),x(-i),x(l1) ,*xs) )/8
    else:
        return   (2*grad(R(1) ,x(-k1),x(-i),x(-j)) 
                 -I(k1,i)*grad(R(1) ,x(l1),x(-l1),x(-j)) 
                 -I(i,j)*grad(R(1)  ,x(-k1),x(l1),x(-l1)) 
                 -I(k1,j)*grad(R(1) ,x(l1),x(-i),x(-l1)) )/8
# hramonics decaying and inner and surfaces 

def h(*idx):
    if idx:
        xs = [x(-i) for i in idx]
        return  grad(R(-1),*xs)
    else:
        return R(-1)
    
def g(*idx):
    if idx:
        return R(len(idx)*2+1) * h(*idx)
    else:
        return 1
def sh(*idx):
    N = -Mul(*[-(3 +2*i) for i in range(len(idx))[:-1]])  
    if idx:
        return R(len(idx)+1) * h(*idx) / N 
    else:
        return 1

    
def n(*idx):
    return x(*idx)*R(-1)

r = symbols("r",positive=True)
# Rto1 = [(R(n),1) for n in range(-20,20)]
# Rto1 = [*Rto1,(r,1)]
# Rtor = [(R(n),r**n) for n in range(-100,100)]
# rtoR = [(r,R(1))]+[(r**n,R(n)) for n in range(-100,100)]
# xtonr = [(x(idx),n(idx)*r) for idx in idxs]

def Rto1(X):
    return X.replace(
        lambda arg: type(arg) == R,lambda arg: 1
        )
def Rtor(X):
    return X.replace(
        lambda arg: type(arg) == R,lambda arg: r**arg.N
        )

def rtoR(X):
    return X.replace(
        lambda arg: arg==r,lambda arg: R(1)
        )
def xtonr(X):
    return X.replace(
            lambda arg: isinstance(arg, Tensor) and arg.head == x,lambda arg: arg * r
        )

def solve_all(Eqs,poly=False):
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
    if poly:
        sols = solve(Poly(TEqs,r).all_coeffs(),const).items()
    else:
        sols = solve(TEqs,const).items()

    for key,s  in sols:
        display(Eq(key, simplify(s)))

    return [(key, factor(s)) for key,s in sols]


def solve_all2(Eqs):
    TEqs=[]
    for eq in Eqs:
        display('Given eq origtens ',eq)
        if isinstance(eq,TensAdd): 
            for e in eq.args:
                if e.args[0] == -1:
                    TEqs.append(e.args[0]*e.args[1])
                else:
                    TEqs.append(e.args[0])
        else:
            TEqs.append(eq.args[0])
        
    TEqsfinal = []
    for eq in TEqs:
        if eq.has(r) :
            expr= factor(eq,r)
            display('Given eq orig ',expr)
            if isinstance(expr,Mul):
                expr= expr.args[-1]
            display('Given eq trans ',expr)
            
            newesq = [expr.coeff(r,ii) for ii in range(-100,100) if expr.coeff(r,ii).has(*const)]
            TEqsfinal += newesq 
        else:
            display('Given eq',eq)
            TEqsfinal.append(eq)


    for e in TEqsfinal:
        display('equations for cst : ',factor(e,r))

    sols = solve(TEqsfinal,const).items()

    for key,s  in sols:
        display(Eq(key,s))
    return [(key, factor(s)) for key,s in sols]



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
Pd3 = Cd[0]*I(k1,k2)+ Cd[1]*g(k1,k2)
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
    if X:
        X =  Rto1(contract_X(X,(I,x,0,0))).doit()
        try: X.canon_bp()
        except:
            return X
        else:
            return X.canon_bp()
    return 0 
        
def compute(X):
    if X:
        X = Rtor(contract_X(X,(I,x,0,0))).doit()
        try: X.canon_bp()
        except:
            return X
        else:
            return X.canon_bp()
    return 0


def compute_Sint(X):
    X = compute_surf(X).int_X((I,x))
    if X: 
        return  simplify(pi * X.doit().contract_delta(I).canon_bp().doit())
    else: 
        return 0 


import subprocess
from sympy.printing.mathematica import mathematica_code
from sympy.parsing.mathematica import parse_mathematica

wolframPath = '/work/fintzin/wolfram/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript'


def max_integrate(F,var,vm,vp,Wolfram=False):
    if Wolfram:
        mathematica_expr_str = mathematica_code(F) 
        minV = mathematica_code(vm)
        maxV = mathematica_code(vp)
        varm = mathematica_code(var)
        wolfram_cmd = (
            f"Integrate[{mathematica_expr_str}, {{{varm},{minV},{maxV}}},"
            +f"Assumptions -> {mathematica_code(phi)} > 0 &&"
            +f"{mathematica_code(lamb)} > 0 &&"
            +f"{mathematica_code(r)} > 0]"
            )
        result = subprocess.run(
                [wolframPath, "-code", wolfram_cmd],
                capture_output=True,
                text=True
            )
        return parse_mathematica(result.stdout.strip())
    else:
        try: integrate(apart(F,var) ,(var,vm,vp))
        except:
            return integrate(F,(var,vm,vp))
        else: 
            return integrate(apart(F,var) ,(var,vm,vp))
            

def compute_Vint(F,vm,vp,Wolfram=False):
    F = simplify(compute(F))
    if isinstance(F,TensAdd):
        terms = []
        for arg in F.args:
            X = sum(1 for n in preorder_traversal(arg)  if getattr(n, 'head', None) == x)
            scalar = Mul(*[a for a in arg.args if not isinstance(a,TensExpr)]) * Mul(*[r]*X)
            tensor = TensMul(*[a for a in arg.args if isinstance(a,TensExpr)]) * R(-X)
            intT = compute_Sint(tensor)
            intS = max_integrate(scalar*r**2,r,vm,vp,Wolfram)
            terms.append(intT*intS)
        return simplify(compute(TensAdd(*terms)))
    
    if isinstance(F,TensMul):
        arg = F
        X = sum(1 for n in preorder_traversal(arg)  if getattr(n, 'head', None) == x)
        scalar = Mul(*[a for a in arg.args if not isinstance(a,TensExpr)]) * Mul(*[r]*X)
        tensor = TensMul(*[a for a in arg.args if isinstance(a,TensExpr)]) * R(-X)
        intT = compute_Sint(tensor)
        intS = max_integrate(scalar*r**2,r,vm,vp,Wolfram)
        return simplify(compute(intT*intS))
    if isinstance(F,Tensor): 
        return  integrate(r**2 ,(r,vm,vp)) * F.int_X((I,x)) *pi

def compute_scalar_int(F,var,vm,vp,Wolfram=False):
    F = simplify(compute(F))
    if isinstance(F,TensAdd):
        terms = []
        for arg in F.args:
            if var ==r:
                X = sum(1 for n in preorder_traversal(arg)  if getattr(n, 'head', None) == x)
                scalar = Mul(*[a for a in arg.args if not isinstance(a,TensExpr)]) * Mul(*[r]*X)
                tensor = TensMul(*[a for a in arg.args if isinstance(a,TensExpr)]) * R(-X)
            else: 
                scalar = Mul(*[a for a in arg.args if not isinstance(a,TensExpr)]) 
                tensor = TensMul(*[a for a in arg.args if isinstance(a,TensExpr)]) 
            intS = max_integrate(scalar,var,vm,vp,Wolfram)
            terms.append(tensor*intS)
        return simplify(compute(TensAdd(*terms)))
    
    if isinstance(F,TensMul):
        arg = F
        if var ==r:
            X = sum(1 for n in preorder_traversal(arg)  if getattr(n, 'head', None) == x)
            scalar = Mul(*[a for a in arg.args if not isinstance(a,TensExpr)]) * Mul(*[r]*X)
            tensor = TensMul(*[a for a in arg.args if isinstance(a,TensExpr)]) * R(-X)
        else: 
            scalar = Mul(*[a for a in arg.args if not isinstance(a,TensExpr)]) 
            tensor = TensMul(*[a for a in arg.args if isinstance(a,TensExpr)]) 
        intS = max_integrate(scalar,var,vm,vp,Wolfram)
        return simplify(compute(tensor*intS))
    if isinstance(F,Tensor): 
        return  integrate(1 ,(var,vm,vp)) * F
    
def series_spe(scalar,X,n,Wolfram):
    if Wolfram:
        mathematica_expr_str = mathematica_code(scalar) 
        wolfram_cmd = (
            f"Normal[Series[{mathematica_expr_str}, {{{mathematica_code(X)},0,{n}}}]]"
            # +f"Assumptions -> {mathematica_code(phi)} > 0 &&"
            # +f"{mathematica_code(lamb)} > 0]"
            )
        result = subprocess.run(
                [wolframPath, "-code", wolfram_cmd],
                capture_output=True,
                text=True
            )
        return parse_mathematica(result.stdout.strip()).subs(X.name,X)
    else: 
        return simplify(series(scalar,X,0,n).removeO()  )

def compute_series(F,X,n=1,Wolfram=False):
    F = compute(F)
    if isinstance(F,TensAdd):
        terms = []
        for arg in F.args:
            number_of_r = sum(1 for n in preorder_traversal(arg)  if getattr(n, 'head', None) == x)
            scalar = Mul(*[a for a in arg.args if not isinstance(a,TensExpr)]) * Mul(*[r]*number_of_r)
            tensor = TensMul(*[a for a in arg.args if isinstance(a,TensExpr)]) / Mul(*[r]*number_of_r)
            series_scalar= series_spe(scalar,X,n,Wolfram)
            terms.append(tensor*series_scalar)
        return simplify(compute(TensAdd(*terms)))
    
    if isinstance(F,TensMul):
        arg = F
        number_of_r = sum(1 for n in preorder_traversal(arg)  if getattr(n, 'head', None) == x)
        scalar = Mul(*[a for a in arg.args if not isinstance(a,TensExpr)]) * Mul(*[r]*number_of_r)
        tensor = TensMul(*[a for a in arg.args if isinstance(a,TensExpr)])  / Mul(*[r]*number_of_r)
        series_scalar= series_spe(scalar,X,n,Wolfram)
        return simplify(compute(tensor*series_scalar))
    
    if isinstance(F,Tensor): 
        return  F
            
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
solsol4 = [
 (Cf[10], Rational(5,48)),
 (Cf[11], Rational(5,48)),
 (Cf[12], -Rational(1,16)),
 (Cf[13], Rational(1,48)),
 (Cf[0], 0),
 (Cf[1], 0),
 (Cf[2], Rational(1,4)),
 (Cf[3], Rational(7,24)),
 (Cf[4], 0),
 (Cf[5], -Rational(1,8)),
 (Cf[6], 0),
 (Cf[7], 0),
 (Cf[8], 0),
 (Cf[9], 0)
]


sol2 = [(Cd[0], -5/(lamb + 1)),
 (Cd[1], -(2*lamb + 3)/(2*(lamb + 1))),
 (Cd[2], -1/(lamb + 1)),
 (Cf[0], (3*lamb + 2)/(2*(lamb + 1))),
 (Cf[1], -(3*lamb + 2)/(4*(lamb + 1))),
 (Cf[2], lamb/(4*(lamb + 1)))]
#solid particles 
solsol2 = [
 (Cf[0], Rational(3,2)),
 (Cf[1], -Rational(3,4)),
 (Cf[2], Rational(1,4))]

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
# solid part
solsol3 = [
 (Cf[0], 0),
 (Cf[1], -Rational(5,3)),
 (Cf[2], Rational(1,2)),
 (Cf[3], -Rational(1,2)),
 (Cf[4], 0),
 (Cf[5], -Rational(1,6))]

sol2B = [(Cd[0], 10/(3*(lamb + 1))),
 (Cd[1], 1/(3*(lamb + 1))),
 (Cd[2], 2/(3*(lamb + 1))),
 (Cf[0], 1/(3*(lamb + 1))),
 (Cf[1], -1/(6*(lamb + 1))),
 (Cf[2], 1/(6*(lamb + 1)))]

sol3B = [
(Cd[0], 0),
(Cd[1], -7/(10*lamb + 10)),
 (Cd[2], (10*lamb + 1)/(60*lamb + 60)),
 (Cd[3], -(10*lamb + 19)/(60*lamb + 60)),
 (Cd[4], 0),#1/(10*lamb + 10)),
 (Cd[5], -1/(12*lamb + 12)),
 (Cf[0],0),
 (Cf[1], -1/(5*lamb + 5)),
 (Cf[2], Rational(1,6)),
 (Cf[3], -Rational(1,6)),
 (Cf[4], 0),
 (Cf[5], -1/(30*lamb + 30))
 ]

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

sol_sol_def2_2 = [
 (Cf[0], Rational( -3,20)),
 (Cf[10], Rational( 1,4)),
 (Cf[11], Rational( 1,4)),
 (Cf[12], Rational( -1,4)),
 (Cf[13], Rational( 1,20)),
 (Cf[1], Rational( -3,20)),
 (Cf[2], Rational( 1,10)),
 (Cf[3], Rational( 1,2)),
 (Cf[4], Rational( 3,40)),
 (Cf[5], Rational( -1,20)),
 (Cf[6], Rational( 3,40)),
 (Cf[7], Rational( -3,20)),
 (Cf[8], Rational( -1,40)),
 (Cf[9], Rational( -1,40))
 ]

sol_sol_def2_3 = [
 (Cf[0], Rational( 16,63)),
 (Cf[10], Rational( 2,5)),
 (Cf[11], Rational( 1,10)),
 (Cf[12], 0),
 (Cf[13], Rational( -1,5)),
 (Cf[14], Rational( -1,10)),
 (Cf[15], Rational( -1,10)),
 (Cf[16], Rational( -2,5)),
 (Cf[17], Rational( 1,5)),
 (Cf[18], Rational( -2,5)),
 (Cf[19], 0),
 (Cf[1], Rational( 1,7)),
 (Cf[20], 0),
 (Cf[21], 0),
 (Cf[22], Rational( 1,10)),
 (Cf[23], Rational( 1,21)),
 (Cf[24], Rational( 1,70)),
 (Cf[25], Rational( 1,70)),
 (Cf[26], Rational( -1,10)),
 (Cf[27], Rational( -3,35)),
 (Cf[28], Rational( -3,35)),
 (Cf[29], Rational( -1,6)),
 (Cf[2], Rational( 1,7)),
 (Cf[30], Rational( 10,63)),
 (Cf[31], Rational( -1,6)),
 (Cf[32], Rational( -1,42)),
 (Cf[3], Rational( 10,21)),
 (Cf[4], Rational( -6,7)),
 (Cf[5], Rational( -6,7)),
 (Cf[6], Rational( -1,3)),
 (Cf[7], Rational( 2,5)),
 (Cf[8], Rational( 1,10)),
 (Cf[9], 0)]



### differentiable scalar functions 

class X(R):
    ### X = r(e.n - 1)/2 
    def __new__(cls,e):
        obj = Symbol.__new__(cls,f'X(\\textbf{{x}},\\textbf{{{e.name}}})')
        obj.e = e
        return obj

    def _eval_partial_derivative(self,v):
        if isinstance(v, TensExpr):
            v=v.substitute_indices(*[(idx,-idx) for idx in v.get_indices()])
            idx = v.get_indices()[0]
            return  (self.e(idx) - v *R(-1))/2
        else:
            return 0
    
    
    def __mul__(self, other):
        return Symbol.__mul__(self, other)

    def __rmul__(self, other):
        return Symbol.__rmul__(self, other)

    def __pow__(self, exponent):
        return Symbol.__pow__(self, exponent)

        
