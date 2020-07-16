#Functions for Normal Shock Waves relations
import numpy as np
import bisection as b
def NSW_TR(g, M):
    TR_NSW = (1 + ((2*g/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2))
    return TR_NSW    #TR_NSW is temp ratio

def NSW_PR(g, M):
    PR_NSW = (1+((2*g)/(g+1))*(M**2-1))
    return PR_NSW    #PR_NSW is pressure ratio

def NSW_ro(g, M):
    ro_NSW = ((g+1)*M**2) /(2+(g-1)*M**2)
    return ro_NSW    #ro_NSW is density ratio

def NSW_M2(g,M):
    M2_NSW = ((1+(((g-1)/2)*(M**2))) / ((g*M**2)-((g-1)/2)))**0.5
    return M2_NSW    #M_2_NSW is M2

def NSW_PO(g,M):
    t = (1 + ((2*g/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2))
    p = (1+((2*g)/(g+1))*(M**2-1))
    dS  = (1005*np.log(t)) - (287*np.log(p))
    DS = (-1*dS)/287
    POr_NSW = np.exp(DS)
    return POr_NSW     #POr_NSW is stagnation pressure ratio

def NSW_DS(g,M):
    t = (1 + ((2*g/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2))   #t is temp ratio used in formula
    p = (1+((2*g)/(g+1))*(M**2-1))    #p is press ratio used in formula
    dS_NSW  = (1005*np.log(t)) - (287*np.log(p))
    return dS_NSW     #dS_NSW is change in entropy
    
def NSW_invPR(Pr_NSW,g):
    M_PR_NSW = (((Pr_NSW-1)*((g+1)/(2*g)))+1)**0.5
    return M_PR_NSW     #M_PR_NSW is Mach no for given press ratio

def NSW_inv_ro(rhor_NSW,g):
    Mro_NSW = ((2*rhor_NSW)/(g+1-rhor_NSW*(g-1)))**0.5
    return Mro_NSW     #Mro_NSW is Mach no for given density ratio

def NSW_invM2(M,g):
    M2_NSW = ((1+(((g-1)/2)*(M**2)))/((g*(M**2))-((g-1)/2)))**0.5
    return M2_NSW     #M_2_NSW is M2

def NSW_invTR(M,x,g):   #x is temp ratio
    M_TR_NSW = ((1 + (((2*g)/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2))-x)
    return  M_TR_NSW      #Mach no for given temp ratio found by calling bisection

def NSW_inv_entropy(M,x,g):   #x is change in entropy
    M_DS_NSW = ((1005*np.log((1 + (((2*g)/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2)))) - (287*np.log(1+((2*g)/(g+1))*(M**2-1))))-x
    return M_DS_NSW         #Mach no for given entropy found by calling bisection

def NSW_invPO(M,x,g):     #x is stag press ratio
    s = (1005*np.log((1 + (((2*g)/(g+1))*(M**2-1)))*((2+((g-1)*M**2))/((g+1)*M**2)))) - (287*np.log(1+((2*g)/(g+1))*(M**2-1)))
    R = (-1*s)/287
    M_POr_NSW = np.exp(R)-x
    return M_POr_NSW          #Mach no for given press ratio found by calling bisection
