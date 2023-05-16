import numpy as np
from .constants import Msun,G,a,c
from .opacity import kappa_lookup
from .config import interp,M,mu,X,Y
from .nuclear import eps_total
from .eos import Pgas,q_eos


def Del_adiabatic(rho,P,T):
    Pg = Pgas(T,rho,mu)
    beta = Pg/P
    ad_inv = (3/2)*beta*(8-7*beta)/(4-3*beta) + 4 - 3*beta
    return 1/ad_inv

def Del_radiative(rho,l,P,T,m):
    kappa = kappa_lookup(T,rho,interp)
    return (3/(16*np.pi*a*c))*P*kappa*l/(T**4*G*m)

def drdxi(r,q):
    return 1/(4*np.pi*r**2*q**2)*(M*Msun)

def dpdxi(r,xi,p):
    return -G*xi/(16*np.pi*r**4*p**3)*(M*Msun)**2

def dFdxi(T,q,X,Y,xi):
    Z = 1-X-Y
    rho = q**2
    epsilon = eps_total(T,rho,X,Y,Z)
    return (epsilon/xi**2)*(M*Msun)

def dTdxi(r,q,xi,F,p,T):    
    m = xi*M*Msun
    rho = q**2
    l = F*xi**2
    P = p**4
    Del_ad = Del_adiabatic(rho,P,T)
    Del_rad = Del_radiative(rho,l,P,T,m)
    if Del_rad <= Del_ad:
        Del = Del_rad
    else:
        Del = Del_ad
    return -G*xi*T/(4*np.pi*r**4*p**4)*Del*(M*Msun)**2

def derivs_xi(ys,xi):
    F,p,r,T = ys
    q = q_eos(p,T,mu)
    dF = dFdxi(T,q,X,Y,xi)
    dp = dpdxi(r,xi,p)
    dr = drdxi(r,q)
    dT = dTdxi(r,q,xi,F,p,T)
    return [dF,dp,dr,dT]

