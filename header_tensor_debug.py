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

def Rto1(X):
    return X.replace(
        lambda arg: isinstance(arg, R),lambda arg: 1
        )
def Rtor(X,r=r):
    return X.replace(
        lambda arg: isinstance(arg, R),lambda arg: r**(arg.N)
        )
def rtoR(X,r=r,x=x):
    return X.replace(
        lambda arg: arg==r,lambda arg: R(1,x=x)
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



def compute_sol(X,sol):
    return simplify(X.canon_bp().subs(sol).doit())

def compute_surf(X):
    if X:
        X =  contract_X(X,(I,x,0,0)).subs(Rto1).doit()
        try: X.canon_bp()
        except:
            return X
        else:
            return X.canon_bp()
    return 0 
        
def compute(X):
    if X:
        X = contract_X(X,(I,x,0,0)).subs(Rtor).doit()
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
        return simplify(series(scalar,x,0,n).removeO()  )

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

        