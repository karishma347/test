
#Functions for isentropic relations

def isnt_TR(g, M):
    T_R= 1/(1+((g-1)/2)*(M**2))
    return T_R   #T_R is temp ratio

def isnt_PR(g,M): 
    P_R = 1/((1+((g-1)/2)*(M**2))**(g/(g-1)))
    return P_R    #P_R is pressure ratio
def isnt_ro(g,M):
    ro_r = 1/((1+((g-1)/2)*(M**2))**(1/(g-1)))
    return ro_r   #ro_r is density ratio

def isnt_AR(g,M):
    A_R = (1/M)*(( (2/(g+1))*(1+((g-1)/2)*M**2))**((g+1)/(2*(g-1))))
    return A_R   #A_R is area ratio

def isnt_invTR(g, T_21):
    M_TR = ((2*((1/T_21)-1))/(g-1))**0.5
    return M_TR  #M_TR is Mach no for given temp ratio

def isnt_invPR(g,P_21):
    M_PR = ((2*(((1/P_21)**((g-1)/g))-1))/(g-1))**0.5
    return M_PR  #M_PR is Mach no for given press ratio

def isnt_inv_ro(g,ro_21):
    M_ro = ((2*(((1/ro_21)**(g-1))-1))/(g-1))**0.5
    return M_ro    #M_ro is Mach no for given density ratio

def isnt_invAR(M,x,g):
    M_Ar = ((1/M)*(( (2/(g+1))*(1+((g-1)/2)*M**2))**((g+1)/(2*(g-1)))))-x
    return M_Ar   #M_Ar is mach no for given area ratio 