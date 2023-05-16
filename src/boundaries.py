import numpy as np
from .constants import Msun,G,a,c,Rsun,sigma_sb
from .config import M,mu,interp,opal_file
from .eos import mu_ionized,q_eos,P_total
from .nuclear import eps_total
from .opacity import interp_opac,kappa_lookup
from scipy.optimize import fminbound

# central boundary conditions

def p_central(M=1,R=1):
    return (1.34e15*(M)**2*(R)**(-4))**(1/4)

def T_central(M,R,mu):
    return 1.15e7*mu*(M)*(R)**(-1)

def p_c(pc_guess,q,xi):
    return (pc_guess**4 - (3*G/8)*((4/3)*np.pi*q**2)**(4/3)*(xi*M*Msun)**(2/3))**(1/4)

def T_c_rad(Tc_guess,q,xi,X,Y,Z):
    rho = q**2
    kappa = kappa_lookup(Tc_guess,rho,interp)
    return (Tc_guess**4 \
            - (1/(2*a*c))*(3/(4*np.pi))**(2/3)*\
            kappa*eps_total(Tc_guess,rho,X,Y,Z)*(q**2)**(4/3)*(xi*M*Msun)**(2/3))**(1/4)

def r_c(q,xi):
    return (3/(4*np.pi*q**2))**(1/3)*(xi*M*Msun)**(1/3)

def F_c(Tc_guess,q,xi,X,Y,Z,psi=1):
    rho = q**2
    return eps_total(Tc_guess,rho,X,Y,Z,psi=1)*M*Msun/xi

def load1(pc_guess,Tc_guess,star_params,xi):
    M,X,Y,Z = star_params
    mu = mu_ionized(X)
    q = q_eos(pc_guess,Tc_guess,mu)
    F = F_c(Tc_guess,q,xi,X,Y,Z,psi=1)
    p = p_c(pc_guess,q,xi)
    r = r_c(q,xi)
    T = T_c_rad(Tc_guess,q,xi,X,Y,Z)
    
    return [F,p,r,T]

# surface boundary conditions

def p_surface(R_guess,kappa):
    return ((2/3)*G*(M*Msun)/((R_guess*Rsun)**2*kappa))**(1/4)

def T_surface(R_guess,F_guess):
    return (F_guess/(4*np.pi*(R_guess*Rsun)**2*sigma_sb))**(1/4)

def delp(pguess,R_guess,Teff,X):
    q = q_eos(pguess,Teff,mu_ionized(X))
    rho = q**2
    kap = kappa_lookup(Teff,rho,interp)
    
    score = np.abs(pguess**4 - p_surface(R_guess,kap)**4)
    return score

def load2(R_guess,F_guess,star_params):
    '''uses pressure bounds from opacity table to find P'''
    M,X,Y,Z = star_params
    Teff = T_surface(R_guess,F_guess)
    _,logRs,logTs,values = interp_opac((X,Y,0,0),opal_file,method="regular",return_values=True)
    Ttabmax = np.power(10,np.amax(logTs))
    Ttabmin = np.power(10,np.amin(logTs))
    rhotabmin = np.power(10,np.amin(logRs)+3*np.amin(logTs)-18)
    rhotabmax = np.power(10,np.amax(logRs)+3*np.amax(logTs)-18)
    P_min = P_total(Ttabmin,rhotabmin,mu_ionized(X))
    p_min = P_min**(1/4)
    P_max = P_total(Ttabmax,rhotabmax,mu_ionized(X))
    p_max = P_max**(1/4)
    F=float(F_guess)
    p=fminbound(delp,x1=p_min,x2=p_max,args=(R_guess,Teff,X))
    r=float(R_guess*Rsun)
    T=float(Teff)
    return [F,p,r,T]

