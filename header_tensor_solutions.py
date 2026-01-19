from  header_tensor_functions import *



###### Solving and physical para

gam, zeta, bet = symbols("Gamma zeta beta")
lamb = symbols('lambda',positive=True)
phi = Symbol('phi', positive=True)

Cf_small = symbols("Cf_{0:50}")
Cd_small = symbols("Cd_{0:50}")
Cf = symbols("Cf_{0:200}")
Cd = symbols("Cd_{0:200}")
Cf2 = symbols("Cf^2_{0:200}")
Cd2 = symbols("Cd^2_{0:200}")
Cf3 = symbols("Cf^3_{0:200}")
Cd3 = symbols("Cd^3_{0:200}")
const_s = [*Cf_small,*Cd_small]

const_all = [*Cf,*Cd,*Cf2,*Cd2,*Cf3,*Cd3]
def const_to_zero(X,const=const_all):
    return X.replace(lambda arg: arg in const, lambda arg: 0 ).doit()

def solve_all(Eqs,cancel_eq=[],const=const_all,Display=True):
    scalar_eqs = []
    for expr in Eqs:
        expr = simplify(expr.canon_bp().subs(cancel_eq).doit())
        if Display: display(expr)
        if isinstance(expr, TensAdd):
            scalar_eqs = scalar_eqs + [Mul(*[a 
                                             for a in arg.args if not isinstance(a,TensExpr)]) 
                                       for arg in expr.args]
        else:
            scalar_eqs.append(Mul(*[a for a in expr.args if not isinstance(a,TensExpr)]))
    

    scalar_eqs_without_r = []
    for eq in scalar_eqs: 
        if eq.has(r): 
            exprs = [a.as_coeff_exponent(r) for a in expand(eq).as_ordered_terms()]
            coeffs = set([a.as_coeff_exponent(r)[1] for a in expand(eq).as_ordered_terms()])
            EQS = [Add(*[c[0] for c in exprs if coeff == c[1]]) for coeff in coeffs]
            EQS = [e for e in EQS if e.has(*const)]
            scalar_eqs_without_r += EQS
        else: 
            scalar_eqs_without_r += [eq]

    if Display: 
        display('system:')
        for eq in scalar_eqs_without_r:
            display(eq)
        display('solutions: ')
    sols = solve(scalar_eqs_without_r,const)
  
    if Display: 
        for key,s  in sols.items(): 
            display(Eq(key ,s) )
    return [(key, factor(s)) for key,s in sols.items()]
##### Iverst the stokes equations for q given force field 

def invert_stokes(forcing,U1=None,P1=None,remove_dirac=True,Display=True,const=const_all):
    ###compute U and P with fourier transform
    if not U1 and not P1:
        FTforcing = compute(FT_expr(compute(forcing)))
        ftU = compute(  FTforcing * (x(-i)*x(k)*R(-4) - Id(-i,k)*R(-2)  ))
        ftP = compute(  FTforcing *  I * x(-i) * R(-2))
        U1 = rtoR(compute(FT_expr(ftU,inverse=True)))
        P1 = rtoR(compute(FT_expr(ftP,inverse=True))) 
        if remove_dirac:
            U1=U1.subs(Dirac(),0).doit()
            P1=P1.subs(Dirac(),0).doit()
        return U1.substitute_indices((k,i)),P1.substitute_indices((k,i))
    else:
        forcingP = grad(rtoR(forcing),x(i))
        Eqp = compute(grad(P1, x(-l1),x(l1)) + forcingP )
        divU = compute(grad(U1,x(i)))
        EqU = compute(grad(U1,x(l1),x(-l1)) - grad(P1,x(-i)) - forcing)
        sols = solve_all([Eqp,divU,EqU],Display=Display,const=const)
        return const_to_zero(U1.subs(sols)),const_to_zero(P1.subs(sols))
    
###### returns the homogeneous solution of a pertubed problem that have particular sol 
def homogeneous_solution(Uho,Pho,Upo,Uhi=None,Phi=None,Upi=None):
    if not (Uhi or Phi or Upi):
        divuf = compute(grad(Uho,x(i)))
        Eqvel2 = compute_surf(Uho + Upo)
        sols = solve_all([divuf,Eqvel2],Display=False)
        Uho = const_to_zero(Uho.subs(sols))
        Pho = const_to_zero(Pho.subs(sols))
        Uho = const_to_zero(Uho.subs(sols))
        Pho = const_to_zero(Pho.subs(sols))
        return Uho,Pho
    divuf = compute(grad(Uho,x(i)))
    divud = compute(grad(Uhi,x(i)))
    Eqvel = compute_surf(Uhi+ Upi - Uho - Upo)
    Eqvel2 = compute_surf((Uhi + Upi)*n(-i) )
    Eqstress = compute_surf(
        n(-j) * E2(Uho + Upo - lamb * (Uhi + Upi)) *(Id(-i,l) - n(-i)*n(l))
    )
    sols = solve_all([divuf,divud,Eqvel,Eqvel2,Eqstress],Display=False)
    Uho = const_to_zero(Uho.subs(sols))
    Pho = const_to_zero(Pho.subs(sols))
    Uhi = const_to_zero(Uhi.subs(sols))
    Phi = const_to_zero(Phi.subs(sols))
    return Uho,Pho,Uhi,Phi

##### Generate harmonics permutations of order n 

def generate_perm_with_cst(ids,Cst):
    i1 = ids[0]
    terms = []
    expr =1 
    for i,ix in enumerate(ids[1:]):
        remaining_id = ids[1:i+1]+ids[i+2:]
        if len(remaining_id) >= 2:
            expr = generate_perm_with_cst(remaining_id,Cst) 
        else:
            expr  = Cst[1][Cst[0]]
            Cst[0] = Cst[0] + 1
        terms.append(Id(i1,ix) * expr)
    return TensAdd(*terms)

from itertools import combinations

def harmonic_perm_I_h(ids,Cst,h0=False,Growing=False):
    odd = False
    if len(ids)%2: odd =True
    terms = []
    if odd: orders = [i+1 for i in range(len(ids))[::2]]
    elif h0: orders  = [0]+ [i+2 for i in range(len(ids))[::2]]
    else: orders = [i+2 for i in range(len(ids))[::2]]
    if Growing: harmo = g
    else: harmo = h

    for order in orders:
        combos = list(combinations(ids, order))
        expr =1 
        for ix in combos:
            remaining_id = [idx for idx in ids if not idx in ix]
            if len(remaining_id) >= 2:
                expr = harmo(*ix)*generate_perm_with_cst(remaining_id,Cst)
            else: 
                expr,Cst[0] = Cst[1][Cst[0]]*harmo(*ix), Cst[0]+1
            terms.append( expr)
    return TensAdd(*terms)

#### not harmonics func 

### non-harmonics function (all permutaioin xII r^n )
def generate_perm_with_r(ids,Cst,Order,logs=False):
    i1 = ids[0]
    terms = []
    expr =1 
    for i,ix in enumerate(ids[1:]):
        remaining_id = ids[1:i+1]+ids[i+2:]
        if len(remaining_id) >= 2:
            expr = generate_perm(remaining_id) 
        r_expr = 0
        for n in range(*Order):
            r_expr += Cst[1][Cst[0]] * R(n)
            Cst[0] = Cst[0]+1
            if logs: 
                r_expr += Cst[1][Cst[0]] * R(n) * log(R(1))
                Cst[0] = Cst[0]+1
        terms.append(Id(i1,ix) * expr * r_expr)
    return TensAdd(*terms)

def perm_I_x(ids,Cst,h0=False,Order=[],logs=False):
    odd = False
    if len(ids)%2: odd =True
    terms = []
    if odd: orders = [i+1 for i in range(len(ids))[::2]]
    elif h0: orders  = [0]+ [i+2 for i in range(len(ids))[::2]]
    else: orders = [i+2 for i in range(len(ids))[::2]]
    for order in orders:
        combos = list(combinations(ids, order))
        expr =1 
        for ix in combos:
            remaining_id = [idx for idx in ids if not idx in ix]
            xis_prod = Mul(*[x(ids) for ids in ix])
            r_expr = 0 
            if len(remaining_id) >= 2:
                expr = xis_prod*generate_perm_with_r(remaining_id,Cst,Order,logs)
            else: 
                for n in range(*Order):
                    r_expr += Cst[1][Cst[0]] * R(n)
                    Cst[0] = Cst[0]+1
                    if logs: 
                        r_expr += Cst[1][Cst[0]] * R(n) * log(R(1))
                        Cst[0] = Cst[0]+1
                expr = xis_prod*r_expr 
            terms.append( expr)
    return TensAdd(*terms)

#### pressure and velocity field harmonics
def P(*idxs,C=[],Growing=False,Order=[],logs=False,h0=False):
    if type(idxs[0]) == int:
        idxs = ik1k2etc[1:idxs[0]]
    Cst = [0,C]
    if not Order:
        expr = harmonic_perm_I_h(idxs,Cst,h0=h0,Growing=Growing)
        return rtoR(Rtor(expr).canon_bp()).doit()
    else: 
        return perm_I_x(idxs,Cst,h0=True,Order=Order,logs=logs)


def U(*idxs,C=[],Growing=False,Order=[],logs=False):
    if type(idxs[0]) == int:
        idxs = ik1k2etc[:idxs[0]]
    Cst = [0,C]
    if not Order: 
        exprp = x(idxs[0])* harmonic_perm_I_h(idxs[1:],Cst,h0=False,Growing=Growing)/2
        expr  = exprp + harmonic_perm_I_h(idxs,Cst,h0=True,Growing=Growing)
        return rtoR(Rtor(expr).canon_bp()).doit()
    else: 
        return perm_I_x(idxs,Cst,h0=True,Order=Order,logs=logs)

# harmonics func 
def P_f(*idxs,Order=[],C=Cf): 
    return P(*idxs,C=C,Order=Order)
def U_f(*idxs,Order=[],C=Cf): 
    return U(*idxs,C=C,Order=Order)
def P_d(*idxs,Order=[],C=Cd,h0=False): 
    return P(*idxs,C=C,Growing=True,Order=Order,h0=h0)
def U_d(*idxs,Order=[],C=Cd): 
    return U(*idxs,C=C,Growing=True,Order=Order)

def E2(Ui):
    return grad(Ui,x(-j))+grad(Ui,x(-j)).substitute_indices((i,j),(j,i))

def Sigma(U,P,mu=1):
    return -P*Id(i,j)+mu*(grad(U,x(-j))+grad(U,x(-j)).substitute_indices((i,j),(j,i)))
def Sigma_f(*idxs):
    return Sigma(U_f(*idxs),P_f(*idxs))




###############################################################
###SOLUTION FOR SINGLE SPHERICAL DROPLET IN BACKGROUND FLOW ###
###############################################################
solD_2 = [
(Cd[0], -5/(lamb + 1)),
 (Cd[1], -(2*lamb + 3)/(2*(lamb + 1))),
 (Cd[2], -1/(lamb + 1)),
 (Cf[0], (3*lamb + 2)/(2*(lamb + 1))),
 (Cf[1], -(3*lamb + 2)/(4*(lamb + 1))),
 (Cf[2], lamb/(4*(lamb + 1)))
 ]
solD_3 = [
 (Cd[0], 7/(2*(lamb + 1))),
 (Cd[2], (2*lamb + 5)/(4*(lamb + 1))),
 (Cd[3], (2*lamb + 5)/(4*(lamb + 1))),
 (Cd[4], 5/(12*(lamb + 1))),
 (Cf[0], -(5*lamb + 2)/(3*(lamb + 1))),
 (Cf[1], 0),
 (Cf[2], 0),
 (Cf[3], 0),
 (Cf[4], -lamb/(6*(lamb + 1)))
] ##cancel stuff 
solD_4 = [
 (Cd[0], (2*lamb - 3)/(2*(lamb + 1))),
 (Cd[10], (lamb - 4)*(2*lamb + 3)/(12*(lamb + 1)*(lamb + 4))),
 (Cd[11], -(2*lamb**2 + 5*lamb + 8)/(12*(lamb + 1)*(lamb + 4))),
 (Cd[12], -(2*lamb**2 + 5*lamb + 8)/(12*(lamb + 1)*(lamb + 4))),
 (Cd[13], -1/(24*(lamb + 1))),
 (Cd[3], -1/(2*(lamb + 1))),
 (Cd[4], -1/(4*(lamb + 1))),
 (Cd[7], (lamb - 1)*(lamb + 5)/(6*(lamb + 1)*(lamb + 4))),
 (Cf[0], lamb/(4*(lamb + 1))),
 (Cf[10], (13*lamb**2 + 10*lamb - 8)/(48*(lamb + 1)*(lamb + 4))),
 (Cf[11], -(lamb - 4)*(3*lamb + 2)/(48*(lamb + 1)*(lamb + 4))),
 (Cf[12], -(lamb - 4)*(3*lamb + 2)/(48*(lamb + 1)*(lamb + 4))),
 (Cf[13], lamb/(48*(lamb + 1))),
 (Cf[3], (7*lamb + 2)/(24*(lamb + 1))),
 (Cf[4], -lamb/(8*(lamb + 1))),
 (Cf[7], (lamb - 1)/(6*(lamb + 1)*(lamb + 4)))
]

###############################################################
### SOLUTION FOR DEFORMING SPHERICAL DROPLET IN STEADY FLOW ###
###############################################################

solD_G_3 = [
 (Cd[0], 7*(2*lamb + 3)/(10*(lamb + 1))),
 (Cd[2], (16*lamb + 19)/(20*(lamb + 1))),
 (Cd[3], (16*lamb + 19)/(20*(lamb + 1))),
 (Cd[4], (2*lamb + 3)/(12*(lamb + 1))),
 (Cf[0], -(19*lamb + 16)/(15*(lamb + 1))),
 (Cf[2], 0),
 (Cf[3], 0),
 (Cf[4], -(3*lamb + 2)/(30*(lamb + 1)))]
 ##this has to be multiplied by  "- grad u_I" 

solD_G_3t2 = [
 (Cd[0], 2*(31*lamb**2 + 28*lamb - 9)/(75*(lamb + 1)**2)),
 (Cd[10], 3*(lamb - 1)*(7*lamb + 8)/(175*(lamb + 1)**2)),
 (Cd[11],-(532*lamb**3 + 266*lamb**2 - 1703*lamb - 1545)/(4200*(lamb + 1)**2)),
 (Cd[12],-(532*lamb**3 + 266*lamb**2 - 1703*lamb - 1545)/(4200*(lamb + 1)**2)),
 (Cd[13], 3*(lamb - 1)*(7*lamb + 8)/(175*(lamb + 1)**2)),
 (Cd[14],-(532*lamb**3 + 3626*lamb**2 + 5647*lamb + 2445)/(4200*(lamb + 1)**2)),
 (Cd[15],-(532*lamb**3 + 3626*lamb**2 + 5647*lamb + 2445)/(4200*(lamb + 1)**2)),
 (Cd[16],(532*lamb**3 + 2870*lamb**2 + 5539*lamb + 3309)/(4200*(lamb + 1)**2)),
 (Cd[17],(532*lamb**3 - 490*lamb**2 - 1811*lamb - 681)/(4200*(lamb + 1)**2)),
 (Cd[19],(532*lamb**3 + 2870*lamb**2 + 5539*lamb + 3309)/(4200*(lamb + 1)**2)),
 (Cd[1], (8*lamb**2 + 119*lamb + 123)/(100*(lamb + 1)**2)),
 (Cd[20],(532*lamb**3 - 490*lamb**2 - 1811*lamb - 681)/(4200*(lamb + 1)**2)),
 (Cd[22],(1708*lamb**3 + 6452*lamb**2 + 4687*lamb - 597)/(9450*(lamb + 1)**2*(2*lamb + 5))),
 (Cd[23], -(56*lamb**2 - 751*lamb - 915)/(7560*(lamb + 1)**2)),
 (Cd[24], -(56*lamb**2 - 751*lamb - 915)/(7560*(lamb + 1)**2)),
 (Cd[25], -(1316*lamb**2 + 2399*lamb + 975)/(7560*(lamb + 1)**2)),
 (Cd[26], -(1316*lamb**2 + 2399*lamb + 975)/(7560*(lamb + 1)**2)),
 (Cd[28],(2*lamb + 3)*(124*lamb + 121)/(2700*(lamb + 1)*(2*lamb + 5))),
 (Cd[29],(2*lamb + 3)*(124*lamb + 121)/(2700*(lamb + 1)*(2*lamb + 5))),
 (Cd[2], (8*lamb**2 + 119*lamb + 123)/(100*(lamb + 1)**2)),
 (Cd[30],7*(2*lamb + 3)*(28*lamb + 97)/(2700*(lamb + 1)*(2*lamb + 5))),
 (Cd[31],7*(2*lamb + 3)*(28*lamb + 97)/(2700*(lamb + 1)*(2*lamb + 5))),
 (Cd[32], 7*(2*lamb + 3)/(2700*(lamb + 1))),
 (Cd[3], -3*(44*lamb**2 + 77*lamb + 29)/(100*(lamb + 1)**2)),
 (Cd[4], -3*(44*lamb**2 + 77*lamb + 29)/(100*(lamb + 1)**2)),
 (Cd[6], 11*(2*lamb + 3)/(270*(lamb + 1))),
 (Cd[8], 3*(lamb - 1)*(7*lamb + 8)/(175*(lamb + 1)**2)),
 (Cd[9], 3*(lamb - 1)*(7*lamb + 8)/(175*(lamb + 1)**2)),
 (Cf[0], 2*(137*lamb**2 + 181*lamb + 32)/(525*(lamb + 1)**2)),
 (Cf[10], 0),
 (Cf[11], -(2*lamb + 3)*(19*lamb + 16)/(300*(lamb + 1))),
 (Cf[12], -(2*lamb + 3)*(19*lamb + 16)/(300*(lamb + 1))),
 (Cf[13], 0),
 (Cf[14], -(2*lamb + 3)*(19*lamb + 16)/(300*(lamb + 1))),
 (Cf[15], -(2*lamb + 3)*(19*lamb + 16)/(300*(lamb + 1))),
 (Cf[16], (2*lamb + 3)*(19*lamb + 16)/(300*(lamb + 1))),
 (Cf[17], (2*lamb + 3)*(19*lamb + 16)/(300*(lamb + 1))),
 (Cf[19], (2*lamb + 3)*(19*lamb + 16)/(300*(lamb + 1))),
 (Cf[1], -2*(269*lamb**2 + 442*lamb + 164)/(525*(lamb + 1)**2)),
 (Cf[20], (2*lamb + 3)*(19*lamb + 16)/(300*(lamb + 1))),
 (Cf[22], (1070*lamb**3 + 4447*lamb**2 + 4619*lamb + 1134)/(4725*(lamb + 1)**2*(2*lamb + 5))),
 (Cf[23], -(410*lamb**2 + 521*lamb + 84)/(4725*(lamb + 1)**2)),
 (Cf[24], -(410*lamb**2 + 521*lamb + 84)/(4725*(lamb + 1)**2)),
 (Cf[25], (125*lamb**2 + 533*lamb + 462)/(9450*(lamb + 1)**2)),
 (Cf[26], (125*lamb**2 + 533*lamb + 462)/(9450*(lamb + 1)**2)),
 (Cf[28],-2*(65*lamb**2 + 251*lamb + 174)/(675*(lamb + 1)*(2*lamb + 5))),
 (Cf[29],-2*(65*lamb**2 + 251*lamb + 174)/(675*(lamb + 1)*(2*lamb + 5))),
 (Cf[2], -2*(269*lamb**2 + 442*lamb + 164)/(525*(lamb + 1)**2)),
 (Cf[30],-7*(2*lamb + 3)*(7*lamb + 4)/(1350*(lamb + 1)*(2*lamb + 5))),
 (Cf[31],-7*(2*lamb + 3)*(7*lamb + 4)/(1350*(lamb + 1)*(2*lamb + 5))),
 (Cf[32], -(143*lamb + 102)/(9450*(lamb + 1))),
 (Cf[3], (127*lamb**2 + 341*lamb + 232)/(525*(lamb + 1)**2)),
 (Cf[4], (127*lamb**2 + 341*lamb + 232)/(525*(lamb + 1)**2)),
 (Cf[6], -(179*lamb + 156)/(675*(lamb + 1))),
 (Cf[8], 0),
 (Cf[9], 0)]
##############################################
###SOLUTION FOR SINGLE SPHERICAL PARTICLES ###
##############################################
sol_2 = [
    (Cf[0], Rational(3,2)), 
    (Cf[1], Rational(-3,4)),
    (Cf[2], Rational(1,4))
]
sol_3 = [
    (Cf[0],Rational( -5,3)), 
    (Cf[2],Rational( -1,2)), 
    (Cf[3],Rational( 1,2)), 
    (Cf[4],Rational( -1,6))
]
sol_4 = [
    (Cf[0],   Rational(1,4)),
    (Cf[10],  Rational(5,48)),
    (Cf[11],  Rational(5,48)),
    (Cf[12],  Rational(-1,16)),
    (Cf[13],  Rational(1,48)),
    (Cf[3],   Rational(7,24)),
    (Cf[4],   Rational(-1,8)),
    (Cf[7],  0)
]

##############################################
### SOLUTION FOR DEFORMED DROPLETS ELLIPSOID #
##############################################

solD_2t2 = [
 (Cd[0], -4*(lamb - 2)/(3*(lamb + 1)**2)),
 (Cd[10], -(23*lamb + 22)/(42*(lamb + 1)*(lamb + 4))),
 (Cd[11], -(23*lamb + 22)/(42*(lamb + 1)*(lamb + 4))),
 (Cd[12], -(lamb + 74)/(21*(lamb + 1)*(lamb + 4))),
 (Cd[13], -1/(21*(lamb + 1))),
 (Cd[1], 2*(lamb - 2)/(lamb + 1)**2),
 (Cd[2], 2*(lamb - 2)/(lamb + 1)**2),
 (Cd[3], -4/(7*(lamb + 1))),
 (Cd[4], -(lamb - 1)/(5*(lamb + 1)**2)),
 (Cd[5], 3*(lamb - 1)/(10*(lamb + 1)**2)),
 (Cd[6], 3*(lamb - 1)/(10*(lamb + 1)**2)),
 (Cd[7], (lamb - 6)*(11*lamb - 61)/(210*(lamb + 1)**2*(lamb + 4))),
 (Cd[8],(197*lamb**2 + 831*lamb - 878)/(420*(lamb + 1)**2*(lamb + 4))),
 (Cd[9],(197*lamb**2 + 831*lamb - 878)/(420*(lamb + 1)**2*(lamb + 4))),
 (Cf[0], (3*lamb**2 - lamb + 8)/(30*(lamb + 1)**2)),
 (Cf[10], (21*lamb**2 + 62*lamb + 52)/(84*(lamb + 1)*(lamb + 4))),
 (Cf[11], (21*lamb**2 + 62*lamb + 52)/(84*(lamb + 1)*(lamb + 4))),
 (Cf[12], -(21*lamb**2 + 22*lamb + 32)/(84*(lamb + 1)*(lamb + 4))),
 (Cf[13], (21*lamb + 4)/(420*(lamb + 1))),
 (Cf[1], -(3*lamb**2 - lamb + 8)/(20*(lamb + 1)**2)),
 (Cf[2], -(3*lamb**2 - lamb + 8)/(20*(lamb + 1)**2)),
 (Cf[3], (7*lamb + 6)/(14*(lamb + 1))),
 (Cf[4], -(3*lamb**2 - lamb + 8)/(60*(lamb + 1)**2)),
 (Cf[5], (3*lamb**2 - lamb + 8)/(40*(lamb + 1)**2)),
 (Cf[6], (3*lamb**2 - lamb + 8)/(40*(lamb + 1)**2)),
 (Cf[7],-(63*lamb**3 + 229*lamb**2 + 438*lamb + 20)/(420*(lamb + 1)**2*(lamb + 4))),
 (Cf[8],-(21*lamb**3 + 143*lamb**2 - 174*lamb + 460)/(840*(lamb + 1)**2*(lamb + 4))),
 (Cf[9],-(21*lamb**3 + 143*lamb**2 - 174*lamb + 460)/(840*(lamb + 1)**2*(lamb + 4)))]

solD_3t2 =[
 (Cd[0], 4*(4*lamb + 1)/(15*(lamb + 1)**2)),
 (Cd[10], 3*(lamb - 1)/(35*(lamb + 1)**2)),
 (Cd[11], (280*lamb**2 + 191*lamb + 19)/(840*(lamb + 1)**2)),
 (Cd[12], (280*lamb**2 + 191*lamb + 19)/(840*(lamb + 1)**2)),
 (Cd[13], 3*(lamb - 1)/(35*(lamb + 1)**2)),
 (Cd[14], -(392*lamb**2 + 1279*lamb + 779)/(840*(lamb + 1)**2)),
 (Cd[15], -(392*lamb**2 + 1279*lamb + 779)/(840*(lamb + 1)**2)),
 (Cd[16], (392*lamb**2 + 1171*lamb + 887)/(840*(lamb + 1)**2)),
 (Cd[17], -(280*lamb**2 + 299*lamb - 89)/(840*(lamb + 1)**2)),
 (Cd[18], (8*lamb + 41)/(28*(lamb + 1)**2)),
 (Cd[19], (392*lamb**2 + 1171*lamb + 887)/(840*(lamb + 1)**2)),
 (Cd[1], (14*lamb**2 + 19*lamb + 17)/(20*(lamb + 1)**2)),
 (Cd[20], -(280*lamb**2 + 299*lamb - 89)/(840*(lamb + 1)**2)),
 (Cd[21], (8*lamb + 41)/(28*(lamb + 1)**2)),
 (Cd[22],(404*lamb**2 + 1495*lamb + 551)/(1890*(lamb + 1)**2*(2*lamb + 5))),
 (Cd[23], (126*lamb**2 + 107*lamb + 89)/(1512*(lamb + 1)**2)),
 (Cd[24], (126*lamb**2 + 107*lamb + 89)/(1512*(lamb + 1)**2)),
 (Cd[25], -(126*lamb**2 + 523*lamb + 289)/(1512*(lamb + 1)**2)),
 (Cd[26], -(126*lamb**2 + 523*lamb + 289)/(1512*(lamb + 1)**2)),
 (Cd[27],-(788*lamb**2 - 1310*lamb - 9523)/(3780*(lamb + 1)**2*(2*lamb + 5))),
 (Cd[28], (124*lamb + 121)/(540*(lamb + 1)*(2*lamb + 5))),
 (Cd[29], (124*lamb + 121)/(540*(lamb + 1)*(2*lamb + 5))),
 (Cd[2], (14*lamb**2 + 19*lamb + 17)/(20*(lamb + 1)**2)),
 (Cd[30], 7*(28*lamb + 97)/(540*(lamb + 1)*(2*lamb + 5))),
 (Cd[31], 7*(28*lamb + 97)/(540*(lamb + 1)*(2*lamb + 5))),
 (Cd[32], 7/(540*(lamb + 1))),
 (Cd[3], -(14*lamb**2 + 51*lamb + 25)/(20*(lamb + 1)**2)),
 (Cd[4], -(14*lamb**2 + 51*lamb + 25)/(20*(lamb + 1)**2)),
 (Cd[5], -(2*lamb - 31)/(6*(lamb + 1)**2)),
 (Cd[6], 11/(54*(lamb + 1))),
 (Cd[7], -(52*lamb + 193)/(210*(lamb + 1)**2)),
 (Cd[8], 3*(lamb - 1)/(35*(lamb + 1)**2)),
 (Cd[9], 3*(lamb - 1)/(35*(lamb + 1)**2)),
 (Cf[0], 2*(25*lamb**2 + 41*lamb + 4)/(105*(lamb + 1)**2)),
 (Cf[10], 0),
 (Cf[11], -(19*lamb + 16)/(60*(lamb + 1))),
 (Cf[12], -(19*lamb + 16)/(60*(lamb + 1))),
 (Cf[13], 0),
 (Cf[14], -(19*lamb + 16)/(60*(lamb + 1))),
 (Cf[15], -(19*lamb + 16)/(60*(lamb + 1))),
 (Cf[16], (19*lamb + 16)/(60*(lamb + 1))),
 (Cf[17], (19*lamb + 16)/(60*(lamb + 1))),
 (Cf[18], 0),
 (Cf[19], (19*lamb + 16)/(60*(lamb + 1))),
 (Cf[1], -2*(52*lamb**2 + 92*lamb + 31)/(105*(lamb + 1)**2)),
 (Cf[20], (19*lamb + 16)/(60*(lamb + 1))),
 (Cf[21], 0),
 (Cf[22],(3*lamb + 4)*(100*lamb**2 + 215*lamb + 7)/(945*(lamb + 1)**2*(2*lamb + 5))),
 (Cf[23], -(81*lamb**2 + 115*lamb + 7)/(945*(lamb + 1)**2)),
 (Cf[24], -(81*lamb**2 + 115*lamb + 7)/(945*(lamb + 1)**2)),
 (Cf[25], (27*lamb**2 + 85*lamb + 112)/(1890*(lamb + 1)**2)),
 (Cf[26], (27*lamb**2 + 85*lamb + 112)/(1890*(lamb + 1)**2)),
 (Cf[27],(180*lamb**3 - 10*lamb**2 - 1241*lamb + 434)/(1890*(lamb + 1)**2*(2*lamb + 5))),
 (Cf[28],-(45*lamb**2 + 110*lamb + 41)/(135*(lamb + 1)*(2*lamb + 5))),
 (Cf[29],-(45*lamb**2 + 110*lamb + 41)/(135*(lamb + 1)*(2*lamb + 5))),
 (Cf[2], -2*(52*lamb**2 + 92*lamb + 31)/(105*(lamb + 1)**2)),
 (Cf[30], -7*(7*lamb + 4)/(270*(lamb + 1)*(2*lamb + 5))),
 (Cf[31], -7*(7*lamb + 4)/(270*(lamb + 1)*(2*lamb + 5))),
 (Cf[32], -(45*lamb + 4)/(1890*(lamb + 1))),
 (Cf[3], (29*lamb**2 + 61*lamb + 50)/(105*(lamb + 1)**2)),
 (Cf[4], (29*lamb**2 + 61*lamb + 50)/(105*(lamb + 1)**2)),
 (Cf[5], (10*lamb**2 - 27*lamb - 4)/(21*(lamb + 1)**2)),
 (Cf[6], -(45*lamb + 22)/(135*(lamb + 1))),
 (Cf[7], 0),
 (Cf[8], 0),
 (Cf[9], 0)
 ]

############################################
### SOLUTION FOR DEFORMED SOLID ELLIPSOID ##
############################################

# XtY means X like the order of the eclt and Y 
# like the order of the deformation which is 
#two for ellipsoids
sol_2t2 = [
    (Cf[0],Rational( 1,10)),
    (Cf[10],Rational( 1,4)),
    (Cf[11],Rational( 1,4)),
    (Cf[12],Rational( -1,4)),
    (Cf[13],Rational( 1,20)),
    (Cf[1],Rational( -3,20)),
    (Cf[2],Rational( -3,20)),
    (Cf[3],Rational( 1,2)),
    (Cf[4],Rational( -1,20)),
    (Cf[5],Rational( 3,40)),
    (Cf[6],Rational( 3,40)),
    (Cf[7],Rational( -3,20)),
    (Cf[8],Rational( -1,40)),
    (Cf[9],Rational( -1,40))
]
sol_3t2 = [
 (Cf[0], Rational( 10,21)),
 (Cf[10], Rational( 1,5)),
 (Cf[11], Rational( -2,5)),
 (Cf[12], Rational( -2,5)),
 (Cf[13], Rational( -1,5)),
 (Cf[14], Rational( -1,10)),
 (Cf[15], Rational( -1,10)),
 (Cf[16], Rational( 2,5)),
 (Cf[17], Rational( 1,10)),
 (Cf[18], Rational( 0)),
 (Cf[19], Rational( 2,5)),
 (Cf[1], Rational( -6,7)),
 (Cf[20], Rational( 1,10)),
 (Cf[21], Rational( 0)),
 (Cf[22], Rational( 10,63)),
 (Cf[23], Rational( -3,35)),
 (Cf[24], Rational( -3,35)),
 (Cf[25], Rational( 1,70)),
 (Cf[26], Rational( 1,70)),
 (Cf[27], Rational( 1,21)),
 (Cf[28], Rational( -1,6)),
 (Cf[29], Rational( -1,6)),
 (Cf[2], Rational( -6,7)),
 (Cf[30], Rational( -1,10)),
 (Cf[31], Rational( 1,10)),
 (Cf[32], Rational( -1,42)),
 (Cf[3], Rational( 1,7)),
 (Cf[4], Rational( 1,7)),
 (Cf[5], Rational( 88,63)),
 (Cf[6], Rational( -1,3)),
 (Cf[7], Rational( 0)),
 (Cf[8], Rational( 0)),
 (Cf[9], Rational( 0))
]
sol_4t2 = [
 (Cf[0],Rational( -3,20)),
 (Cf[100],Rational( -1,112)),
 (Cf[101],Rational( 1,432)),
 (Cf[11],Rational( -1,40)),
 (Cf[14],Rational( -1,40)),
 (Cf[15],Rational( -35,144)),
 (Cf[16],Rational( 55,288)),
 (Cf[17],Rational( 55,288)),
 (Cf[18],Rational( 55,288)),
 (Cf[19],Rational( 55,288)),
 (Cf[1],Rational( 1,4)),
 (Cf[20],Rational( 1,72)),
 (Cf[21],Rational( -5,288)),
 (Cf[22],Rational( -5,288)),
 (Cf[25],Rational( 1,24)),
 (Cf[26],Rational( 3,40)),
 (Cf[27],Rational( -1,8)),
 (Cf[28],Rational( -1,8)),
 (Cf[2],Rational( 1,4)),
 (Cf[30],Rational( 0)),
 (Cf[31],Rational( 0)),
 (Cf[33],Rational( 0)),
 (Cf[34],Rational( 0)),
 (Cf[37],Rational( 1,80)),
 (Cf[40],Rational( 1,80)),
 (Cf[41],Rational( -223,15120)),
 (Cf[42],Rational( 31,3024)),
 (Cf[43],Rational( 31,3024)),
 (Cf[45],Rational( -5,3024)),
 (Cf[46],Rational( -5,3024)),
 (Cf[48],Rational( -5,3024)),
 (Cf[49],Rational( -5,3024)),
 (Cf[4],Rational( 0)),
 (Cf[52],Rational( 7,4320)),
 (Cf[55],Rational( 7,4320)),
 (Cf[56],Rational( -565,6048)),
 (Cf[57],Rational( 1105,12096)),
 (Cf[58],Rational( 1105,12096)),
 (Cf[59],Rational( -565,6048)),
 (Cf[5],Rational( 0)),
 (Cf[60],Rational( 1105,12096)),
 (Cf[61],Rational( 1105,12096)),
 (Cf[62],Rational( 745,12096)),
 (Cf[63],Rational( 745,12096)),
 (Cf[64],Rational( 31,3024)),
 (Cf[65],Rational( 745,12096)),
 (Cf[66],Rational( 745,12096)),
 (Cf[67],Rational( 31,3024)),
 (Cf[68],Rational( 395,6048)),
 (Cf[69],Rational( 205,12096)),
 (Cf[70],Rational( 205,12096)),
 (Cf[71],Rational( -695,12096)),
 (Cf[72],Rational( -155,12096)),
 (Cf[74],Rational( -695,12096)),
 (Cf[75],Rational( -155,12096)),
 (Cf[77],Rational( -695,12096)),
 (Cf[78],Rational( -155,12096)),
 (Cf[7],Rational( 0)),
 (Cf[80],Rational( -695,12096)),
 (Cf[81],Rational( -155,12096)),
 (Cf[83],Rational( -41,3024)),
 (Cf[86],Rational( -1,32)),
 (Cf[87],Rational( 55,4032)),
 (Cf[88],Rational( 55,4032)),
 (Cf[89],Rational( 55,4032)),
 (Cf[8],Rational( 0)),
 (Cf[90],Rational( 55,4032)),
 (Cf[91],Rational( -13,5040)),
 (Cf[92],Rational( -5,4032)),
 (Cf[93],Rational( -5,4032)),
 (Cf[96],Rational( 1,48)),
 (Cf[97],Rational( 1,48)),
 (Cf[98],Rational( 5,336)),
 (Cf[99],Rational( 5,336))
]

##########################################################
### Forced inertial solution around translating sphere ###
##########################################################
u_r =TensorHead('u_r',[L])

U1o = (-1/32)*r**(-3)*1/(lamb**2 + 2*lamb + 1)*(-12*lamb**2*r**3 + 9*lamb**2*r**2 + 3*lamb**2 - 20*lamb*r**3 + 12*lamb*r**2 + 2*lamb - 8*r**3 + 4*r**2)*u_r(i) + (-3/80)*r**(-5)*1/(lamb**2 + 2*lamb + 1)*(15*lamb**2*r**4 - 15*lamb**2 + 25*lamb*r**4 - 13*lamb + 10*r**4 - 2)*u_r(i)*u_r(k)*x(-k) + (-1/320)*r**(-5)*1/(lamb**2 + 2*lamb + 1)*(-90*lamb**2*r**4 + 135*lamb**2*r**3 + 30*lamb**2*r**2 + 15*lamb**2*r - 90*lamb**2 - 150*lamb*r**4 + 180*lamb*r**3 + 38*lamb*r**2 + 10*lamb*r - 78*lamb - 60*r**4 + 60*r**3 + 12*r**2 - 12)*u_r(k)*u_r(-k)*x(i) + (1/32)*r**(-5)*1/(lamb**2 + 2*lamb + 1)*(3*lamb*(3*lamb + 2) - r**2*(9*lamb**2 + 12*lamb + 4))*u_r(k)*x(i)*x(-k) + (3/320)*r**(-7)*1/(lamb**2 + 2*lamb + 1)*(45*lamb**2*r**3 + 15*lamb**2*r - 150*lamb**2 - 10*lamb*r**4*(3*lamb + 5) + 10*lamb*r - 130*lamb - 20*r**4 + r**3*(60*lamb + 20) + r**2*(120*lamb**2 + 158*lamb + 52) - 20)*u_r(k)*u_r(l)*x(i)*x(-k)*x(-l)

U1i = (-1/16)*1/(lamb**2 + 2*lamb + 1)*(3*lamb + r**2*(-6*lamb - 4) + 2)*u_r(i) + (-3/40)*1/(lamb**2 + 2*lamb + 1)*(15*lamb*r**2 - 9*lamb + 10*r**2 - 6)*u_r(i)*u_r(k)*x(-k) + (-1/16)*1/(lamb**2 + 2*lamb + 1)*(3*lamb + 2)*u_r(k)*x(i)*x(-k) + (3/40)*1/(lamb**2 + 2*lamb + 1)*(3*lamb*r**2 - 3*lamb + 2*r**2 - 2)*u_r(k)*u_r(-k)*x(i) + (3/20)*1/(lamb**2 + 2*lamb + 1)*(3*lamb + 2)*u_r(k)*u_r(l)*x(i)*x(-k)*x(-l)

P1i = (-1/16)*r**(-3)*1/(lamb**2 + 2*lamb + 1)*(9*lamb**2 + 12*lamb + 4)*u_r(k)*x(-k) + (-1/160)*r**(-6)*1/(lamb**2 + 2*lamb + 1)*(-150*lamb**2*r**5 + 90*lamb**2*r**4 + 45*lamb**2*r**3 + 30*lamb**2*r**2 + 5*lamb**2 - 250*lamb*r**5 + 120*lamb*r**4 + 83*lamb*r**3 + 20*lamb*r**2 - 100*r**5 + 40*r**4 + 52*r**3)*u_r(k)*u_r(-k) + (1/160)*r**(-8)*1/(lamb**2 + 2*lamb + 1)*(-90*lamb**2*r**5 + 90*lamb**2*r**4 + 135*lamb**2*r**3 + 180*lamb**2*r**2 - 15*lamb**2 - 150*lamb*r**5 + 120*lamb*r**4 + 249*lamb*r**3 + 120*lamb*r**2 - 60*r**5 + 40*r**4 + 156*r**3)*u_r(k)*u_r(l)*x(-k)*x(-l)

P1o = (5/8)*1/(lamb**2 + 2*lamb + 1)*(3*lamb + 2)*u_r(k)*x(-k) + (1/40)*1/lamb*r**2*1/(lamb**2 + 2*lamb + 1)*(63*lamb**2 + 42*lamb + 5*zeta*(r**2 - 1))*u_r(k)*u_r(-k) + (-1/40)*1/lamb*1/(lamb**2 + 2*lamb + 1)*(189*lamb**2 + 126*lamb + 10*r**2*zeta - 15*zeta)*u_r(k)*u_r(l)*x(-k)*x(-l)