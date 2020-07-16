#Bisection for Normal shock waves
import numpy as np
import math
def bisection_NSW(maxmin,funcname,*arglist):   #a and b are guessing values for bisection method
    n_arg_func=funcname.func_code.co_argcount
    if n_arg_func>1:
        arglist_temp0=[maxmin[0]]+list(arglist[0:(n_arg_func-1)]) # creat a list of arguments (pyhton behaviour of counting - exclusive (does not include thelast number given))
        arglist_temp1=[maxmin[1]]+list(arglist[0:(n_arg_func-1)])
    else:
        if n_arg_func==1:
            arglist_temp0=[maxmin[0]]
            arglist_temp1=[maxmin[1]]
#	print(arglist_temp0)
#	print(arglist_temp1)
    arg_array=np.array([[maxmin[0],0],[maxmin[1],0],[0,0]])
    arg_array[0,1]=funcname(*arglist_temp0) # way of passing the values by expansing a list
    arg_array[1,1]=funcname(*arglist_temp1)
#	print(arg_array)
    flag_completion=0
    iterator=1
    tol=1e-10 # This is hard coded will be made optional in the next version
    Max_iter=1e6 # This is hard coded will be made optional in the next version
    while flag_completion==0:
#        
        if (func(a,x,g) * func(b,x,g)>=0):
            print("wrong assumption")
        return
    c=a
    while((b-a)>=0.01):
        c= (a+b)/2
        if(func(c,x,g) == 0.0):
            break
        if(func(c,x,g)*func(a,x,g)<0):
            b=c
        else:
            a=c
    return (a+b)/2

#Bisection for Oblique shock waves
def bisection_OSW(a,b,x,g,func,beta):   #a and b are guessing values for bisection method for OBW relations
    if (func(a,x,g,beta) * func(b,x,g,beta)>=0):
        print("wrong assumption")
        return
    c=a
    while((b-a)>=0.01):
        c= (a+b)/2
        if(func(c,x,g,beta) == 0.0):
            break
        if(func(c,x,g,beta)*func(a,x,g,beta)<0):
            b=c
        else:
            a=c
    return (a+b)/2