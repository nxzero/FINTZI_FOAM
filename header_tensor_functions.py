
###################
## TO DO LIST OPTIMIZE
# we must Rtor at the right place (begining and end not within loops )
from sympy.tensor.tensor import TensorIndexType,tensor_indices, contract_metric,TensorHead, TensAdd, Tensor,TensMul,TensExpr,TensorSymmetry
from sympy.tensor.toperators import PartialDerivative
from sympy import pi
from sympy import *
from sympy.printing.latex import  print_latex

from itertools import product
sym2 = TensorSymmetry.fully_symmetric(2)
asym2 = TensorSymmetry.fully_symmetric(-2)
sym3 = TensorSymmetry.fully_symmetric(3)
sym4 = TensorSymmetry.fully_symmetric(4)

L = TensorIndexType("L",dim=3,dummy_name='L',metric_name='I')
Id = L.delta
d = L.metric
Lev = L.epsilon

x = TensorHead("x", [L]) #position vector
idxs = tensor_indices("i j k l m k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 k_9 l_1 l_2 l_3 l_4 l_5",L)
i, j, k, l, m, k1, k2, k3, k4, k5, k6, k7, k8, k9, l1, l2, l3, l4, l5 = idxs
ik1k2etc = tensor_indices("i k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 k_9 k_10",L)
l1l2l3tc = tensor_indices("i l_1 l_2 l_3 l_4 l_5 l_6 l_7 l_8 l_9 l_10",L)

##fundamental classes of scalar tensor
class Scalar_Tensor(Symbol):
    def __new__(cls):
        obj = Symbol.__new__(cls,f'')
        return obj
    def _eval_partial_derivative(self,v):
        return 0


r = symbols("r")

class R(Scalar_Tensor):
    def __new__(cls,N):
        if N==0:
            return 1
        obj = Symbol.__new__(cls,f'r^{N}')
        obj.N = N
        obj.x = x
        return obj
    
    def _eval_partial_derivative(self,v):
        if isinstance(v, TensExpr):
            v=v.substitute_indices(*[(idx,-idx) for idx in v.get_indices()])
            return (self.N * v * R(self.N-2)).doit()
        else:
            return 0
    
    
class X(Scalar_Tensor):
    ### X = r(e.n - 1)/2 
    def __new__(cls,e,normal_coord=False):
        obj = Symbol.__new__(cls,f'X(\\textbf{{x}},\\textbf{{{e.name}}})')
        obj.e = e
        obj.Re=0
        if normal_coord: 
            obj.Re = Symbol('Re')

        return obj

    def _eval_partial_derivative(self,v):
        if isinstance(v, TensExpr):
            v=v.substitute_indices(*[(idx,-idx) for idx in v.get_indices()])
            idx = v.get_indices()[0]
            if self.Re:
                return  self.Re * (self.e(idx) - v *R(-1))/2
            return (self.e(idx) - v *R(-1))/2
        else:
            return 0


### GRADIENT FUNCTIONs 
def D(expr,v):

    if isinstance(expr, TensExpr):
        if expr.has(Scalar_Tensor):
            if isinstance(expr, TensAdd):
                return TensAdd(*[D(arg, v) for arg in expr.args])
            elif isinstance(expr, TensMul):
                terms = []
                for i, factor in enumerate(expr.args):
                    d_factor = D(factor, v)
                    other_factors = expr.args[:i] + expr.args[i+1:]
                    if d_factor: terms.append(TensMul(d_factor, *other_factors))
                return simplify(TensAdd(*terms))
        else: 
            return expr._eval_partial_derivative(v)
        
    elif isinstance(expr,Scalar_Tensor):
        return expr._eval_partial_derivative(v)
    
    elif isinstance(expr, Add):
        ### this assume that no other scalar are func of the tensor v apart form Scalar_tensor
        return Add(*[D(arg, v) for arg in expr.args if arg.has(Scalar_Tensor)])
    
    elif isinstance(expr, Mul):
        terms = []
        for i, factor in enumerate(expr.args):
            if factor.has(Scalar_Tensor):
                d_factor = D(factor, v)
                other_factors = expr.args[:i] + expr.args[i+1:]
                if d_factor: terms.append(Mul(d_factor, *other_factors))
        return Add(*terms)
    
    elif isinstance(expr,Pow):#derivative of d(f^g)/dx -> (g'ln f + gf'/f) f^g
        return expr.args[0]**expr.args[1]*(
            D(expr.args[1],v)*log(expr.args[0])
            + D(expr.args[0],v)*expr.args[1]/expr.args[0]
        )
    elif isinstance(expr,log):
        return D(expr.args[0],v) / expr.args[0]
    elif isinstance(expr,exp):
        return D(expr.args[0],v) * expr
    elif isinstance(expr, (Integer, Symbol)):
        return S.Zero
    else:
        return S.Zero  # Treat anything else as constant

####the gradient to be used
def grad(expr,*vs):
    for v in vs:
        expr = D(expr, v)
        if isinstance(expr,TensExpr): 
            expr = expr.contract_delta(d)
            if expr: expr = expr.contract_metric(Id) # it takes more time 
            expr = rtoR(simplify(Rtor(expr))) # better for optimisation

    return expr


### Traverse function 
def traverse(expr,func, data=0):
    if isinstance(expr,TensExpr):
        expr = expr.contract_metric(Id)
        if expr: expr = expr.canon_bp() 
        else: return 0
        
        if isinstance(expr, TensAdd):
            return  TensAdd(*[traverse(arg,func, data) for arg in expr.args])
        
        return func(expr,data)
    else: 
        return func(expr,data)

### CONTRACTION FUNCTION 
def contract_X(expr,data):
    X,result = data
    if result == 0: result = R(2)
    if X == 0 or X == x: X = (x,x)

    other_tensors = TensMul(*[a for a in expr.args if not (a.head in X if isinstance(a,TensExpr) else 0)] )
    tensor = [a for a in expr.args if isinstance(a,TensExpr) and a.head in X ] 

    #detect and replace
    contracted_tensor = []
    number_of_contractions = 0 
    for i,ii in enumerate(tensor):
        for j,jj in enumerate(tensor[i+1:]):
            if ii.get_indices()[0] == -jj.get_indices()[0] and ((ii.head, jj.head) == X 
            or (jj.head, ii.head) == X):
                number_of_contractions += 1
                contracted_tensor.append(i)
                contracted_tensor.append(j+i+1)

    final_results = TensMul(*[t for i,t in enumerate(tensor) if i not in contracted_tensor]) * result**number_of_contractions
    return other_tensors*final_results

def contract_expr(expr,data):
    if type(data) == list: 
        for dat in data:
            expr = traverse(expr,contract_X,dat)
    else:
        expr = traverse(expr,contract_X,data) 
    if not expr.get_indices() and isinstance(expr,TensAdd):
        expr = Add(*[a for a in expr.args])
    return expr

###############################""
#### surface and volume integration 
###############################""

def int_X(expr,data):
    ### the goal is xxxx= II + II ... with all the permutation 
    X = data ###the tensor head that represent the position vec
    if expr.has(X):
        tensor = [a for a in expr.args if isinstance(a,TensExpr) and a.head == X ] 
        # * R(len(tensor)), so that tensor represent acctually the normal vector x*R(-1) 
        other_tensors =  R(len(tensor))*TensMul(*[a for a in expr.args 
                                  if not (a.head == X if isinstance(a,TensExpr) else 0)] )
        if len(tensor) % 2: # int of an odd number of n
            return 0
        else:
            factor = Mul(*[i+3 for i in range(len(tensor))[2::2]])*3/4
            exprI = generate_perm(TensMul(*tensor).get_indices())/factor
            return other_tensors*exprI*pi
    else:
        return expr*4*pi ## integration of a cst 

def int_expr(expr,data=x):
    return traverse(expr,int_X,data)


### generate all the permultaions of Id(i,j)*Id(k,l) based on the set of index
def generate_perm(ids):
    i1 = ids[0]
    terms = []
    expr =1 
    for i,ix in enumerate(ids[1:]):
        remaining_id = ids[1:i+1]+ids[i+2:]
        if len(remaining_id) >= 2:
            expr = generate_perm(remaining_id) 
        terms.append(Id(i1,ix) * expr)
    return TensAdd(*terms)


### integrate a scalar 
lamb, E = symbols("lambda E")
phi = Symbol('phi', positive=True)


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
    
#scalar in TensExpre assuming no scalar are hidden in the tensors
def int_scalar(expr,data):
    var,vm,vp,Wolfram = data
    tensors = TensMul(*[a for a in expr.args if isinstance(a,TensExpr)])
    scalars = Mul(*[a for a in expr.args if not isinstance(a,TensExpr)])
    if expr==0: scalars = 0
    return tensors*max_integrate(scalars,*data)

def int_scalar_expr(expr,var,vm=1,vp=oo,Wolfram=False):
    return traverse(expr,int_scalar,(var,vm,vp,Wolfram))


## taylor expansion series
def series_spe(scalar,X,Xo,n,Wolfram):
    if Wolfram:
        mathematica_expr_str = mathematica_code(scalar) 
        wolfram_cmd = (
            f"Normal[Series[{mathematica_expr_str}, {{{mathematica_code(X)},{mathematica_code(Xo)},{n}}}]]"
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
        return series(scalar,X,Xo,n).removeO() 

## series assuming the var is in scalar
def series_tensmul(expr,data):
    if isinstance(expr,TensExpr):
        tensors = TensMul(*[a for a in expr.args if isinstance(a,TensExpr)])
        scalars = Mul(*[a for a in expr.args if not isinstance(a,TensExpr)])
        if expr==0: scalars = 0
        return tensors*series_spe(scalars,*data)
    else:
        return series_spe(expr,*data)

def series_expr(expr,X,Xo=0,n=1,Wolfram=False):
    return traverse(expr,series_tensmul,(X,Xo,n,Wolfram))
## series assuming the var is in scalar
def limit_tensmul(expr,data):
    if isinstance(expr,TensExpr):
        tensors = TensMul(*[a for a in expr.args if isinstance(a,TensExpr)])
        scalars = Mul(*[a for a in expr.args if not isinstance(a,TensExpr)])
        if expr==0: scalars = 0
        return tensors*limit(scalars,*data)
    else: 
        return limit(expr,*data)

def limit_expr(expr,z,z0):
    return traverse(expr,limit_tensmul,(z,z0))

## change of variable r^n --> r^n (1+ n *f* epsion)

def rn_to_frn_tensmul(expr,data):
    X = x
    tensor = [a for a in expr.args if isinstance(a,TensExpr) and a.head == X ] 
    other_tensors = [
        a for a in expr.args if isinstance(a,TensExpr) and not a.head == X 
        ] 
    m = len(tensor)
    scalars = r**m*Mul(*[Rtor(a) for a in expr.args if not isinstance(a,TensExpr)])
    # to apply the change of var r^n = n r^n = r*dr^n /dr   
    scalars = simplify(r*diff(scalars,r)) #maybe apply Rtor all over
    return r**(-m)*TensMul(*(tensor+other_tensors))*scalars

def change_v(expr,f):
    '''changes r--> r*(1+n eps f)'''
    return f *  traverse(expr,rn_to_frn_tensmul,f)

#### ### FOURIER TRANSFORMS AND INVERSE 

class Dirac(Scalar_Tensor):
    ### X = r(e.n - 1)/2 
    def __new__(cls,n=0):
        name = f'\\textbf{{x}}'
        obj = Symbol.__new__(cls,f'\delta({name})')
        if n:
            obj = Symbol.__new__(cls,f'\\nabla^{({n})}\delta({name})')
        obj.order = n 
        obj.X = X
        obj.Nabla = TensorHead(f'\\nabla',[L])
        return obj

    def _eval_partial_derivative(self,v):
        if isinstance(v, TensExpr):
            idxs = [-idx for idx in v.get_indices()]
            return TensMul(*[self.Nabla(idx) for idx in idxs])*self
        else:
            return 0

def FT_expr(expr,inverse=False,X=x,r_scalar=r,r_func=R):
    if expr.has(r_func): expr = Rtor(expr)
    return traverse(expr,FT_tensor,(inverse,X,r_scalar,r_func))

def FT_tensor(expr,data):
    inverse,X,r_scalar,r_func = data 
    if isinstance(expr,TensExpr):
        tensor = TensMul(*[
            a for a in expr.args if isinstance(a,TensExpr) and not (a.head == X or a.head == Dirac().Nabla) 
            ]) 
        tensor_x = [a for a in expr.args if isinstance(a,TensExpr) and a.head == X ] 
        tensor_nabla = [a for a in expr.args if isinstance(a,TensExpr) and a.head == Dirac().Nabla] 
        scalars = Mul(*[a for a in expr.args if not isinstance(a,TensExpr)])
        if expr==0: scalars = 0   

        # scalars = apart(scalars,r_scalar)
        ### transform scalars 
        FTscalar = rtoR(FT_scalar(scalars,data)) #need to be diff

        ### transform tensor just the nabla that comes with the dirac
        FTtensor_nabla = TensMul(
            *[X(n.get_indices()[0]) for n in tensor_nabla]
            ) * I**len(tensor_nabla)
        if inverse: FTtensor_nabla =FTtensor_nabla*(-1)**len(tensor_nabla)
        
        ## goes from x to nabla
        n = len(tensor_x)
        xs = [X(-x.get_indices()[0]) for x in tensor_x ]
        FTexpr = I**n*grad(FTscalar*FTtensor_nabla,*xs)

        if inverse: FTexpr = FTexpr * (-1)**n 

        FTexpr = FTexpr * tensor 
        

        if isinstance(FTexpr,TensExpr):
            return compute(FTexpr)
        else:
            return FTexpr
    else: 
        return FT_scalar(expr,data)

def FT_scalar(expr,data):
    inverse,X,r_scalar,r_func = data
    FTexpr = 0
    for term in expand(expr).as_ordered_terms():
        if not term.has(Dirac) and not term.has(ln(r_scalar)):
            coef,expo = term.as_coeff_exponent(r_scalar) 
            FTexpr += coef * FT(expo,data)
        elif term.has(Dirac):
            coef = Mul(*[a for a in term.args if not isinstance(a,Dirac) ])
            dirdel = Mul(*[a for a in term.args if isinstance(a,Dirac) ])
            if isinstance(term,Dirac): dirdel = term
            FTexpr += coef* FT(dirdel,data)
        elif term.has(ln(r_scalar)):
            coef = Mul(*[a for a in term.args if not a.has(r_scalar) ])
            dirdel = Mul(*[a for a in term.args if a.has(r_scalar)])
            if isinstance(term,ln): dirdel = term
            FTexpr += coef* FT(dirdel,data)
    return simplify(FTexpr)

def FT(n,data): #returns the fourier transform ---or inverse--- of r^n
    inverse,X,r_scalar,r_func = data
    if (n%2==0 and n<0 )or (n%2==1 and n>0):
        FT_r_to_n = (2**(3+n) * pi**Rational(3,2)
                     *gamma(Rational(3,2)+Rational(n,2))
                     *r_scalar**(-n-3))/gamma(-Rational(n,2))
        
    elif (n%2==1 and n<0 and n != -1):
        h = (-n-3)/2
        FT_r_to_n = 4*pi*(-1)**h *r_scalar**(2*h)*(
            -log(r_scalar)
            +log(2)
            +digamma(h+Rational(3,2))/2
            +digamma(h+1)/2
        )/gamma(2+2*h)

    elif n==(-1): 
        FT_r_to_n = 4*pi*r_scalar**(-2) 

    elif n%2==0 and n>=0:
        h = n/2
        FT_r_to_n = 8*pi**3*(-1)**(h) * Dirac(n)
    
    elif isinstance(n,Dirac):  
        n = n.order
        FT_r_to_n = I**(n) * r_scalar**(n) 

    elif n.has(ln(r_scalar)):
        if isinstance(n,ln):  
            FT_r_to_n = - 2 * pi**2 * R(-3) +8*pi**3* Dirac() *(1-EulerGamma)
        else: 
            if isinstance(n,Pow): return display('POW')
            rs = Mul(*[a for a in n.args if not isinstance(a,ln)])
            lnr = Mul(*[a for a in n.args if isinstance(a,ln)])
            coef,n = rs.as_coeff_exponent(r_scalar) 
            if lnr == ln(r_scalar): 
                if (n%2==0 and n<0 ) or (n%2==1 and n>0):#r^{-2}ln r r^{-4}ln r
                    FT_r_to_n = simplify(
                        4*2**n*pi**(Rational(3,2))*r_scalar**(-n - 3)*(
                        -2*ln(r) + polygamma(0, -Rational(n,2)) 
                        + polygamma(0, Rational(n,2) + Rational(3,2)) + ln(4)
                        )*gamma(Rational(n,2) + Rational(3,2))/gamma(-Rational(n,2))
                    )
                elif (n%2==1 and n<0 and n != -1): #r^{-3}ln r r^{-5}ln r
                    n_var=Symbol('n_var')
                    h = (-n_var-3)/2
                    FT_r_to_n = -(-1)**(-n/2 + Rational(-3, 2))*pi*r_scalar**(-n - 3)*(
                        -4*log(r_scalar)**2 
                        + 2*log(r_scalar)*polygamma(0, -n/2) 
                        + 4*log(r_scalar)*polygamma(0, -n - 1) 
                        + 2*log(r_scalar)*polygamma(0, -n/2 + Rational(-1, 2)) 
                        + 4*log(2)*log(r_scalar) 
                        - 2*I*pi*log(r_scalar) 
                        - 2*polygamma(0, -n/2)*polygamma(0, -n - 1) 
                        + I*pi*polygamma(0, -n/2) 
                        - 2*polygamma(0, -n - 1)*polygamma(0, -n/2 + Rational(-1, 2)) 
                        - 4*log(2)*polygamma(0, -n - 1) 
                        + I*pi*polygamma(0, -n/2 + Rational(-1, 2)) 
                        + polygamma(1, -n/2) 
                        + polygamma(1, -n/2 + Rational(-1, 2)) 
                        + 2*I*pi*log(2))/gamma(-n - 1)
                elif n==(-1): 
                    FT_r_to_n = 4*pi*r_scalar**(-2)*(-ln(r_scalar)-EulerGamma) 
                elif n%2==0 and n>=0: # r^2 ln (r) r^4 ln(r) and so on 
                    FT_r_to_n = (- 2*pi**2* gamma(2+n)*r_scalar**(-n-3)/(-1)**Rational(n,2)
                                 +8*pi**3 * (-1)**Rational(n,2) * Dirac(n) *(
                                     ln(2) +polygamma(0,Rational(n,2) +Rational(3,2))/2
                                     +polygamma(0,Rational(n,2) +1)/2
                                 ))

    if inverse: 
        return FT_r_to_n/(8*pi**3)
    else: 
        return FT_r_to_n



###FUNDAMENTAL HARMONICS

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
    
def n(idx):
    return x(idx)*R(-1)
### simplify function 

def Rtor(X):
    if X: return X.replace(lambda arg: type(arg) == R,lambda arg: r**arg.N)
    else: return 0
def Rto1(X):
    return Rtor(X).subs(r,1)
def Rto2(X,at_r):
    return Rto1( xtonR(Rtor(X)).subs(r,at_r) )

def rtoR(X):
    return X.replace(lambda arg:  arg==r, lambda arg: R(1)
            ).replace(lambda arg:  type(arg)==Pow and type(arg.args[0]) == R, 
            lambda arg: R(arg.args[1]))
def xtonR(X):
    return X.replace(lambda arg:  isinstance(arg,Tensor) and arg.head==x, 
                     lambda arg: arg*R(-1)*r)

####alll the simplifying functions 
def simplify_tens(X):
    return simplify(X.canon_bp().contract_delta(Id))

def compute_sol(X,sol):
    return simplify(X.subs(sol))

def compute_surf(X,at_r = 1):
    if X:
        return  simplify(Rto2(contract_expr(X,(x,0)),at_r=at_r))
    return 0 
        
def compute(X, contract_rules= (x,0),scalars_rules=Rtor):
    if X:
        return simplify(scalars_rules(contract_expr(X,contract_rules)))
    return 0

##integration over the unit surface 
def S_int(X):
    return  simplify(compute_surf(int_expr(X,x)))

## itegration on a volume around the particle
def V_int(expr,x=x,vm=1,vp=oo,Wolfram=False): ## the var is on r,x 
    # the r**2 is the metric coef
    surf_int = Rtor(int_expr(expr,x))
    # if surf_int: #if it is not zero
    return simplify(int_scalar_expr(r**2 *surf_int,r,vm,vp,Wolfram))
    # else: 
        # return 0


## integration over a def vol dV = (1+ 3 f )dV^0
# it return only the perturbed part 
def V_intf(expr,f,x=x,vm=1,vp=oo,Wolfram=False):
    return V_int(change_v(expr,f) + expr*3*f, 
                 x=x,
                 vm=vm,
                 vp=vp,
                 Wolfram=Wolfram)

def S_intf(expr,f):
    return S_int(change_v(expr,f) +  expr*2*f)