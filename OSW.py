#Functions for Oblique Shock waves Relations
import math
def OSW_Mn1(M1, g, beta):
    b = beta*math.pi/180
    Mn_1_OSW = M1*math.sin(b)
    return Mn_1_OSW    #Mn_1_OSW is Mn1

def OSW_Mn2(M1,g,beta):
    Mn_1 = OSW_Mn1(M1, g, beta)
    Mn_2_OSW = ((1+(((g-1)/2)*(Mn_1**2))) / ((g*Mn_1**2)-((g-1)/2)))**0.5
    return Mn_2_OSW    #Mn_2_OSW is Mn2

def OSW_M2(M1,g,beta, theta):
    Mn2 = OSW_Mn2(M1,g,beta)
    b = beta*math.pi/180
    t = theta*math.pi/180
    M2_OSW = Mn2/math.sin(b-t)
    return M2_OSW    #M2_OSW is M2 

def OSW_PR(M1, g,beta):
    Mn_1 = OSW_Mn1(M1,g,beta)
    PR_OSW = (1+((2*g)/(g+1))*(Mn_1**2-1))
    return PR_OSW    #PR_OSW is pressure ratio

def OSW_ro(M1, g,beta):
    Mn_1 = OSW_Mn1(M1,g,beta)
    ro_OSW = ((g+1)*Mn_1**2) /(2+(g-1)*Mn_1**2)
    return ro_OSW    #ro_OSW is density ratio

def OSW_TR(M1,g,beta):
    Mn_1 = OSW_Mn1(M1,g,beta)
    TR_OSW = (1 + ((2*g/(g+1))*(Mn_1**2-1)))*((2+((g-1)*Mn_1**2))/((g+1)*Mn_1**2))
    return TR_OSW    #TR_OSW is temp ratio

def OSW_invPR(Pr_OSW,g):
    M1n_PR_OSW = (((Pr_OSW-1)*((g+1)/(2*g)))+1)**0.5
    return M1n_PR_OSW     #M1n_PR_OSW is Mach no for given press ratio

def OSW_inv_ro(rhor_OSW,g):
    M1n_ro_OSW = ((2*rhor_OSW)/(g+1-rhor_OSW*(g-1)))**0.5
    return M1n_ro_OSW     #M1n_ro_OSW is Mach no for given density ratio

def OSW_invM2(M1n,g):
    M2n_OSW = ((1+(((g-1)/2)*(M1n**2)))/((g*(M1n**2))-((g-1)/2)))**0.5
    return M2n_OSW     #M2n_OSW is M2n

def OSW_invTR(M1,x,g,beta):   #x is temp ratio using bisection
    Mn_1 = OSW_Mn1(M1,g,beta)
    M1n_TR_OSW = ((1 + (((2*g)/(g+1))*(Mn_1**2-1)))*((2+((g-1)*Mn_1**2))/((g+1)*Mn_1**2))-x)
    return M1n_TR_OSW    #mach no is found by calling bisection
