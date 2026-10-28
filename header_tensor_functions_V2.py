###################
## TO DO LIST OPTIMIZE
# For instance if we pass just a tensor argument 
# to all func it doesn't work but that should not happen 
# Write all things with replace with wild expression wild tensor ? 
# print functions 
import time
import sys 

def printt( str, *str2 ): 
    print(time.ctime(),'|' ,   str, *str2   )
    sys.stdout.flush()
    
from sympy.tensor.tensor import TensorIndexType,tensor_indices, contract_metric,TensorHead, TensAdd, Tensor,TensMul,TensExpr,TensorSymmetry
from itertools import combinations, permutations
from sympy import *
import pickle 
import re 
# metric 
L = TensorIndexType("L",dim=3,dummy_name='L',metric_name='I')
Id = L.delta
d = L.metric
Lev = L.epsilon

x = TensorHead("x", [L]) #position vector
nt = TensorHead("n", [L]) #normal vector

# symbols("r",positive=True) #norm of position with respect to origin 


i, j, k, l, m, o,p = tensor_indices("i j k l m o p",L)
k1,k2,k3,k4,k5,k6,k7,k8,k9 = tensor_indices("k_1:10",L)
l1,l2,l3,l4,l5,l6,l7,l8,l9 = tensor_indices("l_1:10",L)

##fundamental classes of scalar tensor such as r 
class Scalar_Tensor(Symbol):
    def __new__(cls):
        obj = Symbol.__new__(cls,f'')
        return obj
    def _eval_partial_derivative(self,v):
        return 0



class R(Scalar_Tensor):
    def __new__(cls):
        obj = Symbol.__new__(cls,f'R',positive=True)
        return obj
    
    def _eval_partial_derivative(self,v):# here v represent x since r = (x.x)^(1/2)
        if isinstance(v, TensExpr):
            v=v.substitute_indices(*[(idx,-idx) for idx in v.get_indices()])
            return  v  / R() 
        else:
            return 0
# we will be using only the instace "r" of the class "R"    
r = R()

class X(Scalar_Tensor):
    ### X = r(e.n - 1)/2 
    def __new__(cls,e,normal_coord=False):
        obj = Symbol.__new__(cls,f'X')
        obj.e = e
        obj.Re = 0
        if normal_coord: 
            obj.Re = Symbol('Re',positive=True)

        return obj

    def _eval_partial_derivative(self,v):
        if isinstance(v, TensExpr):
            v=v.substitute_indices(*[(idx,-idx) for idx in v.get_indices()])
            idx = v.get_indices()[0]
            if self.Re:
                return  self.Re * (self.e(idx) - v/r )/2
            return (self.e(idx) - v/r )/2
        else:
            return 0
### simplify function with r

def Rto1(X,at_r=1):
    return X.subs(r,at_r)

def xton(X):
    return expand(X.replace(lambda arg:  isinstance(arg,Tensor) and arg.head.name=='x', 
                     lambda arg: r*nt(*arg.get_indices())))
def ntox(X):
    return X.replace(lambda arg:  isinstance(arg,Tensor) and arg.head.name=='n', 
                     lambda arg: x(*arg.get_indices())/r ).doit()

def x_to_mx(X):
    '''this is usefull to get the opposite if F(x) i.e. F(-x) '''
    if isinstance(X,list): return [x_to_mx(X) for X in X]
    if X:
        return X.replace(
            lambda arg: isinstance(arg,Tensor) and arg.head.name == 'x',
            lambda arg: -arg
        ).doit()
    else: 
        return X

# tensor that may posses a derivative idicated with :
# \grad_j T^ijk = T^ijk_j (it just add indexes but that is enough)
# the symmetry is like it is symmetric on the derivative and that T has a 
# non modified sym
class TensorFieldHead(TensorHead):
    def __call__(self, *args):    
        # this must be re-coded using head.symmetry.base and head.rank 
        return TensorField(self, args)
    
class TensorField(Tensor):
    def _eval_partial_derivative(self, v):
        idx = -v.get_indices()[0]
        sym = self.head.sym + [self.head.rank - self.head.order +1 ] 
        TF = TensorFieldHead(self.head.name +"'",
                             self.head.index_types+[L],
                             TensorSymmetry.direct_product(*sym)
                             )(*self.get_indices(),idx)
        TF.head.sym = self.head.sym
        TF.head.order = self.head.order
        return TF 

### DERIVATIVE ON TENSOR EXPR FUNCTIONs 
def tensor_derivative(expr,v):
    if isinstance(expr, TensExpr):
        if expr.has(Scalar_Tensor):
            if isinstance(expr, TensAdd):
                return TensAdd(*[tensor_derivative(arg, v) for arg in expr.args])
            elif isinstance(expr, TensMul):
                terms = []
                for i, factor in enumerate(expr.args):
                    d_factor = tensor_derivative(factor, v)
                    other_factors = expr.args[:i] + expr.args[i+1:]
                    if d_factor: terms.append(TensMul(d_factor, *other_factors))
                return expand(TensAdd(*terms))
        else: 
            return expr._eval_partial_derivative(v)
        
    elif isinstance(expr,Scalar_Tensor):
        return expr._eval_partial_derivative(v)
    
    elif isinstance(expr, Add):
        ### this assume that no other scalar are func of the tensor v apart form Scalar_tensor
        return Add(*[tensor_derivative(arg, v) for arg in expr.args if arg.has(Scalar_Tensor)])
    
    elif isinstance(expr, Mul):
        terms = []
        for i, factor in enumerate(expr.args):
            if factor.has(Scalar_Tensor):
                d_factor = tensor_derivative(factor, v)
                other_factors = expr.args[:i] + expr.args[i+1:]
                if d_factor: terms.append(Mul(d_factor, *other_factors))
        return Add(*terms)
    
    elif isinstance(expr,Pow):#derivative of d(f^g)/dx -> (g'ln f + gf'/f) f^g
        return expr.args[0]**expr.args[1]*(
            tensor_derivative(expr.args[1],v)*log(expr.args[0])
            + tensor_derivative(expr.args[0],v)*expr.args[1]/expr.args[0]
        )
    elif isinstance(expr,log):
        return tensor_derivative(expr.args[0],v) / expr.args[0]
    elif isinstance(expr,exp):
        return tensor_derivative(expr.args[0],v) * expr
    elif isinstance(expr, (Integer, Symbol)):
        return S.Zero
    else:
        return S.Zero  # Treat anything else as constant

####the gradient to be used
def grad(expr,*vs):
    for v in vs:
        expr = tensor_derivative(expr, v) 
        if isinstance(expr,TensExpr): 
            expr = expr.contract_delta(d)
            if expr: expr = expr.contract_metric(Id) # it takes more time 
            expr = expand(expr) # better for optimisation
    return expr
 

### CONTRACTION FUNCTION 
def contract_TensMul(expr,X,result):
    other_tensors = TensMul(*[
        a for a in expr.args 
        if not (a.head in X if isinstance(a,TensExpr) else 0)] )
    tensor = [a for a in expr.args if isinstance(a,TensExpr) and a.head in X ] 
    #detect and replace
    contracted_tensor = []
    number_of_contractions = 0 
    #this code works for contraction with only one idx 
    for i,ii in enumerate(tensor):
        for j,jj in enumerate(tensor[i+1:]):
            if ii.head.rank == 1 and jj.head.rank == 1:
                if (ii.get_indices()[0] == -jj.get_indices()[0] 
                    and ((ii.head, jj.head) == X or (jj.head, ii.head) == X)):
                    number_of_contractions += 1
                    contracted_tensor.append(i)
                    contracted_tensor.append(j+i+1)
            if ii.head.rank == 2 or jj.head.rank == 2: 
                if ((-jj.get_indices()[0] in ii.get_indices() 
                      or -ii.get_indices()[0] in jj.get_indices())
                    and ((ii.head, jj.head) == X or (jj.head, ii.head) == X)
                    and result==0):
                    number_of_contractions += 1
                    contracted_tensor.append(i)
                    contracted_tensor.append(j+i+1) 

    final_results = TensMul(*[
        t for i,t in enumerate(tensor) if i not in contracted_tensor
        ]) * result**number_of_contractions
    expr = other_tensors*final_results
    if isinstance(expr,TensExpr):
        if not expr.get_indices():
            expr = Mul(*[a for a in expr.args])
    return expr
    
def contract_expr(expr,data):    
    if not isinstance(expr,TensExpr): return expr 
    if isinstance(data,list): 
        for dat in data: expr = contract_expr(expr,dat)
        return expr  
    
    X,result = data

    expr =  expr.contract_metric(Id).replace(
        lambda arg: isinstance(arg,TensMul) 
        and arg.has(X[0]) and arg.has(X[1])
        and arg._get_dummy_indices_set(),   
        lambda arg: contract_TensMul(arg,X,result)
    )
    
    if isinstance(expr,TensExpr):
        if not expr.get_indices() and isinstance(expr,TensAdd):
            expr = Add(*[a for a in expr.args])
    return expr

###############################""
#### Integration on the unit sphere (on the angles ) 
###############################""

def int_expr(expr):
    """This function iter through every terms of an expression and integrate it 
    on the unit sphere surface \iint f(\theta,\phi) d\theta d\phi """
    if isinstance(expr,Tensor) and expr.head == x: 
        return 0
    elif isinstance(expr,Tensor):
        return 4*pi*expr 
    expr = expr.contract_metric(Id).replace(
        lambda arg: isinstance(arg,TensMul) and arg.count(x) == 0,   
        lambda arg: 4*pi*arg
    ).replace(
        lambda arg: isinstance(arg,TensMul) and arg.count(x) != 0,   
        lambda arg: int_X(arg) if (not arg.count(x) % 2 and not arg.has(TensorField)) 
        else arg if arg.has(TensorField) else 0 
     )
    return expr

def generate_perm_for_Id(ids):
    """this recursive funcitons generate all the combianaison of index representing the 
    indices of a product of identity tensor"""
    terms,expr,i1 = [],[],ids[0]
    for i,ix in enumerate(ids[1:]):
        remaining_id = ids[1:i+1]+ids[i+2:]
        if len(remaining_id) >= 2:
            expr = generate_perm_for_Id(remaining_id) 
            if len(expr) > 1: terms += [[(i1,ix)] +  exp for exp in expr ]
            else: terms += [[(i1,ix)] + expr]    
        else: terms += [(i1,ix)] + expr   
    return terms 

def get_tensor_symmetry(tensor): 
    """That function retuns the groups of symmetry of a tensor but 
    does't specify if it is symmetric or skew symmetric"""
    if isinstance(tensor,TensMul): 
        syms = []
        for ten in tensor.args: 
            if get_tensor_symmetry(ten):
                syms = syms + get_tensor_symmetry(ten) 
        return syms 
    
    syms = tensor.head.symmetry.base
    ids = tensor.get_indices()
    ids_sym = [] 
    for ii in range(len(syms)): 
        if syms[ii]  == (syms[ii-1]+1): 
            ids_sym = ids_sym[:-1] + [(*ids_sym[-1],ids[syms[ii]+1])]
        else:
            ids_sym = ids_sym + [(ids[syms[ii]],ids[syms[ii]+1])]
    for ii in ids: 
        if ii not in [jj for pair in ids_sym for jj in pair]: 
            ids_sym = ids_sym + [(ii,)]
    return ids_sym

def int_X(expr):
    """This function replace the Tensor *xxxx --> the unique permutaions of Tensor *Id*Id 
    This is optimized in the same way as it is done in the REDUCE program of :
    On Computer Generated Analytic Solutions to ... (bathèse-biesel et al 1972)"""
    idx_xs = TensMul(*[a for a in expr.args 
                       if isinstance(a,TensExpr) and a.head == x]).get_free_indices()     
    tensor = TensMul(*[a for a in expr.args if isinstance(a,TensExpr) and not a.head == x])
    if len(tensor.args)==1: tensor = tensor.args[0]
    factor = Mul(*[i+3 for i in range(len(idx_xs))[2::2]])*3/4/pi
    scalars = r**(len(idx_xs))*Mul(*[a for a in expr.args if not isinstance(a,TensExpr)])/factor
    if len(idx_xs) == 2: 
            return scalars*tensor*Id(*idx_xs)
        
    # get the list of contracted index and others 
    idx_contracted = set([-i for i in tensor.get_free_indices()]) & set(idx_xs)
    idx_not_contracted = tuple(set(idx_xs) - idx_contracted)
    idx_syms =  [tuple(-i for i in turp if -i in idx_contracted) 
                 for turp in get_tensor_symmetry(tensor)]
    idx_syms = idx_syms + [(i,) for i in idx_not_contracted]  
    
    # label the symmetric index by number and generate unique permutations, and get occurence of each
    scheme_idx_xs = [[ii for ii in range(len(idx_syms)) if idx in idx_syms[ii]][0] for idx in idx_xs]
    perms = [tuple(sorted(perm)) for perm in generate_perm_for_Id(scheme_idx_xs)]
    setperms = set(perms)
    reduced_perms = [[[i for i  in pair] for pair in perm] for perm in setperms]
    occurences = [perms.count(perm) for perm in setperms]

    # replace the markes by the actual index 
    for ii in range(len(reduced_perms)):
        markers = [0]*len(idx_syms)
        for jj in range(len(reduced_perms[ii])): 
            for idx in range(len(reduced_perms[ii][jj])):
                rank = reduced_perms[ii][jj][idx]
                reduced_perms[ii][jj][idx] = idx_syms[rank][markers[rank]]
                markers[rank] += 1 
                
    # built the sum, contract, and symplify 
    return TensAdd(*[
        occ*scalars*(tensor*TensMul(*[Id(i,j)for i,j in perm])).contract_delta(Id).canon_bp() 
        for occ, perm in zip(occurences,reduced_perms)]).doit() 
    
### generate all the permultaions of Id(i,j)*Id(k,l) based on the set of index
def generate_perm_Id(ids):
    if len(ids) == 2: return Id(*ids)
    else: return TensAdd(*[TensMul(*[
        Id(*pair) for pair in perm]) for perm in generate_perm_for_Id(ids)])
    
#########################################################################
# Below functions to generate explicitely spherical harmonics (a lot faster 
# than doing explicit differentiation )
#########################################################################

def generate_perm_for_Id_x(ids):
    '''This gives the list of all permutation of the type x^(m) * Id^(n) 
    These fundamentales permutation are used to compute explicitely the spherical harmonics 
    of the form grad^(n) 1/r or even grad^(n) r'''
    if len(ids)%2: orders = [i+1 for i in range(len(ids))[::2]]
    else : orders  = [0]+ [i+2 for i in range(len(ids))[::2]]
    terms = []
    for order in orders[::-1]:
        combos = list(combinations(ids, order))
        for ix in combos:
            remaining_id = ids.copy()
            for idx in ix: remaining_id.remove(idx)
            if len(remaining_id) >= 2:
                if len(remaining_id) >= 4:
                    expr = [[ix]  + deltas for deltas in  generate_perm_for_Id(remaining_id) ]
                    for exp in expr: terms.append(exp)
                else:  terms.append([ix] + generate_perm_for_Id(remaining_id) )
            else:  terms.append([ix])
    return terms                  

def hr(expr, order_r,*idxs): 
    if expr: 
        tensor = TensMul(*[a for a in expr.args if isinstance(a,TensExpr) and not a.head == x])
        if len(tensor.args)==1: tensor = tensor.args[0]
        scalars = Mul(*[a for a in expr.args if not isinstance(a,TensExpr)])
        # get the list of contracted index and others 
        idx_contracted = set([-i for i in tensor.get_free_indices()]) & set(idxs)
        idx_not_contracted = tuple(set(idxs) - idx_contracted)
        idx_syms =  [tuple(-i for i in turp if -i in idx_contracted)  
                     for turp in get_tensor_symmetry(tensor)]
        idx_syms = idx_syms + [(i,) for i in idx_not_contracted]  
    else: 
        scalars,tensor =1,1
        scheme_idx_xs = list(idxs)
        idx_syms = []+ [(i,) for i in idxs]  
    
    # label the symmetric index by number and generate unique permutations, and get occurence of each
    scheme_idx_xs = [[ii for ii in range(len(idx_syms)) if idx in idx_syms[ii]][0] for idx in idxs]
    perms = [ tuple([perm[0]] + sorted(perm[1:]))  for perm in generate_perm_for_Id_x(scheme_idx_xs)]
    setperms = set(perms)
    reduced_perms = [[[i for i  in pair] for pair in perm] for perm in setperms]
    occurences = [perms.count(perm) for perm in setperms]
    
    # replace the markes by the actual index 
    for ii in range(len(reduced_perms)):
        markers = [0]*len(idx_syms)
        for jj in range(len(reduced_perms[ii])): 
            for idx in range(len(reduced_perms[ii][jj])):
                rank = reduced_perms[ii][jj][idx]
                reduced_perms[ii][jj][idx] = idx_syms[rank][markers[rank]]
                markers[rank] += 1 

    # writing the final results with the correct coeff 
    results = []
    for perm,occ in zip(reduced_perms , occurences): 
        res = (scalars *tensor * TensMul(
                *[x(i) for i in perm[0]],
                *[Id(i,j) for i,j in perm[1:]]
                )).contract_delta(Id).canon_bp()
        order = res.count(x) + (len(idxs)-res.count(x))//2
        # order_r cannot be other thing than -1 and 1 because this coef is not yet generalized 
        coef_2 = - 1 if order_r == 1 else 0 if order_r == -1 else 0.666
        coef = occ * facto_facto(order+ coef_2 )* r**(-2*order+order_r)  
        results += coef * res 
    return results 

#################################################################
### Below is a pattern to apply scalar function on the terms ###
### just in front of the tensors, this includes limit, integration (scalar), 
### and other. It also contain mathematica 
#################################################################
def scalar_operation_TensMul(expr,scalar_func,*data):
    tensors = TensMul(*[a for a in expr.args if isinstance(a,TensExpr)])
    scalars = Mul(*[a for a in expr.args if not isinstance(a,TensExpr)])
    if expr==0: scalars = S.Zero
    return tensors*scalar_func(scalars,*data)

def scalar_operation(expr,scalar_func,*data):
    if isinstance(expr,Tensor): return expr*scalar_func(1,*data)
    if not isinstance(expr,TensExpr): return scalar_func(expr, *data)
    return expr.canon_bp().replace(
        lambda arg: isinstance(arg,TensMul),
        lambda arg: scalar_operation_TensMul(arg,scalar_func,*data),
    )    
    
lamb, E = symbols("lambda E")
phi = Symbol('phi', positive=True)


import subprocess
from sympy.printing.mathematica import mathematica_code
from sympy.parsing.mathematica import parse_mathematica

wolframPath = '/work/fintzin/wolfram/SystemFiles/Kernel/Binaries/Linux-x86-64/wolframscript'

# scalar integration using mathematica 
def max_integrate(F,var,vm,vp,Wolfram):
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

## taylor expansion on sclar including mathematica
def max_series(scalar,X,Xo,n,Wolfram):
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
## below is short-cut for these funciton and limits 

def int_scalar_expr(expr,var,vm=1,vp=oo,Wolfram=False):
    return scalar_operation(expr, max_integrate,var,vm,vp,Wolfram)
def series_expr(expr,X,Xo=0,n=2,Wolfram=False):
    return scalar_operation(expr, max_series,X,Xo,n,Wolfram)
def limit_expr(expr,z,z0):
    return scalar_operation(expr, limit,z,z0)

def approxR(X,n,Xo=oo):
    '''Since r = (x . x)^1/2 its series expansion isn't quite teh same as scalar so here is a 
    short cut to take series expr on r variable  
    By defaults this takes the taylor expansion with respect to the distance function R at oo
    for decaying harmonics but it could also be aplly near 0 '''
    if n < 0 and Xo == oo: n = 0
    if isinstance(X,TensExpr):
        return expand(ntox(series_expr(xton(X),r,Xo,n)))
    else:
        return expand(series_expr(X,r,Xo,n))
    
## change of variable r^n --> r^n (1+ n *f* epsion)
def change_v(expr,f):
    '''changes r--> r*(1+n eps f) and so r^n -> 1 + r^n n eps +exp^2 '''
    # for that we use n r^n = r *diff(r^n,n)
    return f * r* ntox(scalar_operation(xton(expr),diff,r))

######################################################################
#### FOURIER TRANSOFMR AND INVERSE FOURIER TRANS
#######################################################""""""""

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

def FT(n,inverse): #returns the fourier transform ---or inverse--- of r^n
    if (n%2==0 and n<0 )or (n%2==1 and n>0):
        FT_r_to_n = (2**(3+n) * pi**Rational(3,2)
                     *gamma(Rational(3,2)+Rational(n,2))
                     *r**(-n-3))/gamma(-Rational(n,2))
        
    elif (n%2==1 and n<0 and n != -1):
        h = (-n-3)/2
        FT_r_to_n = 4*pi*(-1)**h *r**(2*h)*(
            -log(r)
            +log(2)
            +digamma(h+Rational(3,2))/2
            +digamma(h+1)/2
        )/gamma(2+2*h)

    elif n==(-1): 
        FT_r_to_n = 4*pi*r**(-2) 

    elif n%2==0 and n>=0:
        h = n/2
        FT_r_to_n = 8*pi**3*(-1)**(h) * Dirac(n)
    
    elif isinstance(n,Dirac):  
        n = n.order
        FT_r_to_n = I**(n) * r**(n) 

    elif n.has(ln(r)):
        if isinstance(n,ln):  
            FT_r_to_n = - 2 * pi**2 * r**(-3) +8*pi**3* Dirac() *(1-EulerGamma)
        else: 
            if isinstance(n,Pow): return display('POW')
            rs = Mul(*[a for a in n.args if not isinstance(a,ln)])
            lnr = Mul(*[a for a in n.args if isinstance(a,ln)])
            coef,n = rs.as_coeff_exponent(r) 
            if lnr == ln(r): 
                if (n%2==0 and n<0 ) or (n%2==1 and n>0):#r^{-2}ln r r^{-4}ln r
                    FT_r_to_n = expand(
                        4*2**n*pi**(Rational(3,2))*r**(-n - 3)*(
                        -2*ln(r) + polygamma(0, -Rational(n,2)) 
                        + polygamma(0, Rational(n,2) + Rational(3,2)) + ln(4)
                        )*gamma(Rational(n,2) + Rational(3,2))/gamma(-Rational(n,2))
                    )
                elif (n%2==1 and n<0 and n != -1): #r^{-3}ln r r^{-5}ln r
                    n_var=Symbol('n_var')
                    h = (-n_var-3)/2
                    FT_r_to_n = -(-1)**(-n/2 + Rational(-3, 2))*pi*r**(-n - 3)*(
                        -4*log(r)**2 
                        + 2*log(r)*polygamma(0, -n/2) 
                        + 4*log(r)*polygamma(0, -n - 1) 
                        + 2*log(r)*polygamma(0, -n/2 + Rational(-1, 2)) 
                        + 4*log(2)*log(r) 
                        - 2*I*pi*log(r) 
                        - 2*polygamma(0, -n/2)*polygamma(0, -n - 1) 
                        + I*pi*polygamma(0, -n/2) 
                        - 2*polygamma(0, -n - 1)*polygamma(0, -n/2 + Rational(-1, 2)) 
                        - 4*log(2)*polygamma(0, -n - 1) 
                        + I*pi*polygamma(0, -n/2 + Rational(-1, 2)) 
                        + polygamma(1, -n/2) 
                        + polygamma(1, -n/2 + Rational(-1, 2)) 
                        + 2*I*pi*log(2))/gamma(-n - 1)
                elif n==(-1): 
                    FT_r_to_n = 4*pi*r**(-2)*(-ln(r)-EulerGamma) 
                elif n%2==0 and n>=0: # r^2 ln (r) r^4 ln(r) and so on 
                    FT_r_to_n = (- 2*pi**2* gamma(2+n)*r**(-n-3)/(-1)**Rational(n,2)
                                 +8*pi**3 * (-1)**Rational(n,2) * Dirac(n) *(
                                     ln(2) +polygamma(0,Rational(n,2) +Rational(3,2))/2
                                     +polygamma(0,Rational(n,2) +1)/2
                                 ))

    if inverse: 
        return FT_r_to_n/(8*pi**3)
    else: 
        return FT_r_to_n
 
# TODO add the fourier transform of exp and its inverse ! 
def FT_scalar(expr,inverse):
    FTexpr = 0
    for term in expand(expr).as_ordered_terms():
        if not term.has(Dirac) and not term.has(ln(r)):
            coef,expo = term.as_coeff_exponent(r) 
            FTexpr += coef * FT(expo,inverse)
        elif term.has(Dirac):
            coef = Mul(*[a for a in term.args if not isinstance(a,Dirac) ])
            dirdel = Mul(*[a for a in term.args if isinstance(a,Dirac) ])
            if isinstance(term,Dirac): dirdel = term
            FTexpr += coef* FT(dirdel,inverse)
        elif term.has(ln(r)):
            coef = Mul(*[a for a in term.args if not a.has(r) ])
            dirdel = Mul(*[a for a in term.args if a.has(r)])
            if isinstance(term,ln): dirdel = term
            FTexpr += coef* FT(dirdel,inverse)
    return expand(FTexpr)


def FT_TensMul(expr,inverse):
    tensor = TensMul(*[
        a for a in expr.args if isinstance(a,TensExpr) 
        and not (a.head == x or a.head == Dirac().Nabla) 
        ]) 
    tensor_x = [a for a in expr.args if isinstance(a,TensExpr) and a.head == x ] 
    tensor_nabla = [a for a in expr.args if isinstance(a,TensExpr) and a.head == Dirac().Nabla] 
    scalars = Mul(*[a for a in expr.args if not isinstance(a,TensExpr)])
    if expr==0: scalars = S.Zero   
    ### transform scalars 
    FTscalar = FT_scalar(scalars,inverse) #need to be diff
    ### transform tensor just the nabla that comes with the dirac
    FTtensor_nabla = TensMul(
        *[x(n.get_indices()[0]) for n in tensor_nabla]
        ) * I**len(tensor_nabla)
    if inverse: FTtensor_nabla =FTtensor_nabla*(-1)**len(tensor_nabla)
    
    ## goes from x to nabla
    n = len(tensor_x)
    xs = [x(-xx.get_indices()[0]) for xx in tensor_x ]
    FTexpr = I**n*grad(FTscalar*FTtensor_nabla,*xs)
    if inverse: FTexpr = FTexpr * (-1)**n 
    FTexpr = FTexpr * tensor 
    
    if isinstance(FTexpr,TensExpr):
        return compute(FTexpr)
    else:
        return FTexpr
    
def FT_expr(expr,inverse=False):
    '''fourier transform an expression based on that x is the position vector and r its norm'''
    if isinstance(expr,Tensor) and not expr.head == x: return expr * FT_scalar(1,inverse)
    if isinstance(expr,Tensor) and expr.head == x: return  -FT_expr(-expr,inverse) # the two is completely arbitrary there
    if not isinstance(expr,TensExpr): return FT_scalar(expr,inverse)
    return compute(expr).replace(
        lambda arg: isinstance(arg,TensMul),
        lambda arg: expr * FT_scalar(1,inverse)
        if not expr.has(x,r) else FT_TensMul(arg,inverse)
    )

####################################""
### FUNDAMENTAL HARMONICS

def h(*idx):
    '''decaying harmonics'''
    return hr(None , -1, *idx)
    
def g(*idx):
    '''growing harmonics'''
    return r**(len(idx)*2+1) * h(*idx)
    
def sh(*idx):
    '''surface harmonics'''
    N = -Mul(*[-(3 +2*i) for i in range(len(idx))[:-1]])  
    if idx:
        return r**(len(idx)+1) * h(*idx) / N 
    else:
        return 1

def n(idx):
    return x(idx)/r # norm of position 

def facto_facto(n):
    return (-1)**n * Mul(*[(2*n -2*k-1) for k in range(n)])

#############################################
### Standart Simplifying functions 

def compute_surf(X,contract_rules= ((x,x),r**2),at_r = 1):
    if X:
        return  expand(Rto1(contract_expr(X,contract_rules),at_r=at_r))
    return 0 
        
def compute(X, contract_rules= ((x,x),r**2)):
    if X:
        return expand(contract_expr(X,contract_rules))
    return 0


#############################################
###INTEGRATION OVER THE UNIT SPHERE

##integration over the unit surface 
def S_int(X,at_r = 1): #it acctually rpz r^2 dS 
    return  expand(compute_surf(r**2*int_expr(X),at_r=at_r))

## itegration on a volume around the particle
def V_int(expr ,vm=0,vp=1,Wolfram=False): ## the var is on r,x 
    # the r**2 is the metric coef
    surf_int = r**2 * int_expr(expr)
    # if surf_int: #if it is not zero
    return expand(int_scalar_expr(surf_int,r,vm,vp,Wolfram))
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

## integration over a deformed surface dS = (1+ 2 f )dS^0
# it return only the perturbed part  prop to f
def S_intf(expr,f):
    return S_int(change_v(expr,f) +  expr*2*f)


###########################################
# Substitutions of tensors expression 

def replace_tesnor(Xs , names , Ts):
    '''replace a Tensor by an expr of same order. 
    This becomes quite complicated when replacing a tensor contracted with others
    this is why teh substitution is complexe here'''
    if isinstance(Xs,list): 
        return [replace_tesnor(X,names,Ts) for X in Xs]
    
    if not isinstance(Ts,list): Ts = [Ts]
    
    if isinstance(Xs,TensExpr): 
        return compute(Xs.replace(
            lambda arg: isinstance(arg,Tensor) and arg.head.name in names, 
            lambda arg: Ts[names.index(arg.head.name)].substitute_indices(
                    *zip(sorted(Ts[names.index(arg.head.name)].get_free_indices() )
                    if isinstance(Ts[names.index(arg.head.name)].get_free_indices(),list) 
                    else [Ts[names.index(arg.head.name)].get_free_indices()],
                    arg.get_indices()
                    if isinstance(arg.get_indices(),list)
                    else [arg.get_indices()]),
                *zip(tensor_indices(f"L_0:100",L), tensor_indices(f"l_0:100",L))
                ) if isinstance(Ts[names.index(arg.head.name)],TensExpr) 
            else Ts[names.index(arg.head.name)]
        ))
    else:
        return Xs
        

def convert_to_scalar(X):
    '''This get tensor to scalar symbols.'''
    if isinstance(X,list): return [convert_to_scalar(X) for X in X]
    X=X.replace(
        lambda arg: isinstance(arg,Tensor),
        lambda arg: symbols(arg.head.name+str(''.join([str(i) for i in arg.get_indices()])))
    )
    if isinstance(X,TensAdd):
        return Add(*[Mul(*mul.args) for mul in X.args ])
    if isinstance(X,TensMul):
        return Mul(*X.args) 


# going in the other way
def symbol_to_tensor(X,parameters_names):
    '''Given tensors name, and order and symmetry we can reconstruct
    expression based on scalar expression generated with convert to scalar funciton
    for instace it is very case specific 
    
    names are the list of parameters (scalars that we know are tensor in solutions)'''
    if X.name[0] in parameters_names:
        tokens = re.findall(r'-?[A-Za-z](?:_\d+)?', X.name )[1:]
        # the first token is the name of the tensor not idx 
        idxsstr = " ".join(t.lstrip('-').replace('L_', 'l') for t in tokens)
        idxs = tensor_indices(idxsstr,L)
        if not isinstance(idxs,list): idxs = [idxs]
        idxs = [(a if tokens[i][0] != '-' else -a ) for i,a in enumerate(idxs) ]
        o = len(idxs)
        Gn = TensorHead(f'{X.name[0]}_{o}',[L]*o,
                        TensorSymmetry.direct_product(1,o-1))
        return Gn(*idxs)
    
    if X.name[0] == 'd':
        tokens = re.findall(r'-?[A-Za-z](?:_\d+)?', X.name )[1:]
        idxsstr = " ".join(t.lstrip('-').replace('L_', 'l')  for t in tokens)
        idxs = tensor_indices(idxsstr,L)
        if not isinstance(idxs,list): idxs = [idxs]
        idxs = [(a if tokens[i][0] != '-' else -a ) for i,a in enumerate(idxs) ]
        d = TensorHead('d',[L])
        return d(*idxs)
    
    if X.name[0] == 'x':
        tokens = re.findall(r'-?[A-Za-z](?:_\d+)?', X.name )[1:]
        idxsstr = " ".join(t.lstrip('-').replace('L_', 'l')  for t in tokens)
        idxs = tensor_indices(idxsstr,L)
        if not isinstance(idxs,list): idxs = [idxs]
        idxs = [(a if tokens[i][0] != '-' else -a ) for i,a in enumerate(idxs) ]
        return x(*idxs)
    
    if X.name[0] == 'V':
        name = re.findall(r'-?[A-Za-z](?:_\d+)?', X.name )[0]
        tokens = re.findall(r'-?[A-Za-z](?:_\d+)?', X.name )[1:]
        idxsstr = " ".join(t.lstrip('-').replace('L_', 'l')  for t in tokens)
        idxs = tensor_indices(idxsstr,L)
        if not isinstance(idxs,list): idxs = [idxs]
        idxs = [(a if tokens[i][0] != '-' else -a ) for i,a in enumerate(idxs) ]
        return TensorHead(name,[L]*len(idxs))(*idxs)
    
    if X.name[:2] == 'KD':
        tokens = re.findall(r'-?[A-Za-z](?:_\d+)?', X.name[2:])
        idxsstr = " ".join(t.lstrip('-').replace('L_', 'l') for t in tokens)
        idxs = tensor_indices(idxsstr,L)
        return Id(*idxs)

def scalar_to_tensor(X,parameters_names):
    return X.replace(
        lambda arg: isinstance(arg,Symbol),
        lambda arg: symbol_to_tensor(arg,parameters_names)
        )

    
##########################################################################
# BELOW SOME SPECIFIQUE FUNCION FOR TENSOR EQUAITON SOLUTIONS ############
# These are originally defined in DDB_methd.ipy 

def Ltraces(expr,names=['P','U']):
    '''take the trace to zero on any index of the tensor that posses the names in the
    list '''
    if isinstance(expr,list): return [Ltraces(X) for X in expr]
    if isinstance(expr,TensExpr):
        return expr.contract_delta(Id).replace(# remove tensor for divergence free
            lambda arg: (isinstance(arg,Tensor) and arg.head.name[0] in names
                        and sorted(list(set( [ii if ii.is_up else -ii for ii in arg.get_indices()] )))
                        != sorted([ii if ii.is_up else -ii for ii in arg.get_indices()])),
            lambda arg: 0
        ).doit()
    else: return expr

def GFtraces(expr,names=['G','M']):#
    '''These tensor represent for example grad^(n) U such that any contraction of
    the first idx with any other is zero so that div(U) = 0.
    Also this include the fact that grad^2 grad^{(n)} G = \grad^{n+1}1/r and therefore
    in the case of contracted this becomes fully symmetric.
    And finally it is zero for Lap^2 grad^n G = 0.
    The properties of GF can be reflected into the tensors that are contracted with it
    and so additional names can be inclued'''
    if isinstance(expr,list): return [GFtraces(X,names) for X in expr]
    if isinstance(expr,TensExpr):
        return expr.contract_delta(Id).canon_bp().replace(# remove tensor for divergence free
        lambda arg: (isinstance(arg,Tensor) and arg.head.name[0] in names
                     and - arg.get_indices()[0] in arg.get_indices()[1:]),
        lambda arg: 0
        ).doit().replace(# remove any double contracted tensor
        lambda arg: (isinstance(arg,Tensor) and arg.head.name[0] in names
                     and len(arg.get_indices())-len(arg.get_free_indices()) > 3),
        lambda arg: 0
        ).doit().replace(
        # When laplacian is actived the tensor becomes fully sym because grad grad p is symetric
        # however we have to make sure that the dummy idx does not get contracted with the velocity vector
        # do not apply this on div u terms --> they are canceled at the step before
        lambda arg: (isinstance(arg,Tensor) and arg.head.name[0] in names and arg.head.rank > 3 
                     and len(arg.get_indices()) - len(arg.get_free_indices()) == 2), # G may be remove
        lambda arg: 
            TensorHead(arg.head.name,
                        [L]*arg.head.rank,
                        TensorSymmetry.direct_product(len(arg.get_free_indices()),len(arg._get_dummy_indices_set()))
                        )(*arg.get_free_indices(),*arg._get_dummy_indices_set())
        ).canon_bp() #the doit() are NEEDED !!!
    else: return expr
    
    
def remove_traces(expr,GF_names = ['G','M'],L_names=['P','X','F','U']): 
    return GFtraces(Ltraces(expr, names=L_names), names=GF_names)


def GF_permutes(expr,names=['G','M']):
    '''Generates permutaion of \lapacian grad u etc to explicit the symetry constrains
    This is usefull for exmple when substituing \Lap \grad U by a tensor (which may not be explicit 
    so that the result becomes explicit). 
    
    WHY THIS FUNCTION ISN'T YET NEEDED FOR STOKES FLOW FIEDLDS: 
    if the expression is defined by symbolic tensor such as G or V and its grad this funciton 
    would just generate the bi-laplacian terms in addition to the div u terms. Both vanish 
    either ways..... 
    '''
    if isinstance(expr,TensExpr):
        return expr.contract_delta(Id).canon_bp().replace(# when laplacian is actived the rest is symmetric because \grad^2 G = \grad\grad 1/r or something
        lambda arg: (isinstance(arg,Tensor) and arg.head.name[0] in names and arg.head.rank > 3 
                     and len(arg.get_indices()) - len(arg.get_free_indices()) == 2), # G may be remove
        lambda arg: TensAdd(*[TensorHead(arg.head.name,arg.head.rank*[L])(*arg.get_indices())] 
                            + [TensorHead(arg.head.name,arg.head.rank*[L])(*arg.get_indices()).substitute_indices(*idx) 
                              for idx in 
                        [(pair, tuple(reversed(pair))) 
                        for pair in combinations(arg.get_free_indices(),2) 
                        if arg.get_free_indices()[0] in pair ]]
                        ) / len(arg.get_free_indices())
        ).canon_bp()
    else:
        return expr

def removeKD(X,idxs = i):     
    if isinstance(X,list): return [removeKD(smallX,idxs) for smallX in X]
    if isinstance(idxs,list): 
        for idx in idxs: X = removeKD(X,idx)
        return X
    return X.replace(
        lambda arg: isinstance(arg,Tensor) and arg.head == Id 
        and (idxs in arg.get_indices()),
        lambda arg: 0
    ).replace(
        lambda arg: isinstance(arg,TensMul) and arg.count(Id) >= 2,
        lambda arg: 0  
    ).doit()


def sym(X):
    return (X/2+X.substitute_indices((i,j),(j,i))/2).canon_bp()