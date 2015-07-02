"""Katharin Jensen
This code implements the 4th order Runge-Kutta method for the dynamic analysis equations"""

#imports
import math


def RK4(u0,u1,y0,v0,H,L,k,ea,m,g):
    #inputs: u at the beginning (u0) and end (u1) of the step
    #        y and v at the beginning of the time step
    #        parameters H, L, k, ea, and m
    #outputs: y and v at the end of the time step
    
    #find step size h
    h = u1-u0
    
    #find next y value
    K_1y = v0
    K_2y = v0 + (h/2)*K_1y
    K_3y = v0 + (h/2)*K_2y
    K_4y = v0 + h*K_3y
    
    y1 = y0 + (h/6)*(K_1y + 2*K_2y + 2*K_3y + K_4y)
    
    #find next v value
    #ensure it doesn't get math errors from too large steps
    if ((u0)-(y0))/(H-L) <= -1:
        K_v1 = -(k/(4*m))*(L-(y0)) + (ea/m)*math.atanh(-.9999999999999999) + g
    elif ((u0)-(y0))/(H-L) >= 1:
        K_v1 = -(k/(4*m))*(L-(y0)) + (ea/m)*math.atanh(.9999999999999999) + g
    else:
        K_v1 = -(k/(4*m))*(L-(y0)) + (ea/m)*math.atanh(((u0)-(y0))/(H-L)) + g
    
    if ((u0+(h/2))-(y0+(h/2)*K_v1))/(H-L) <= -1:
        K_v2 = -(k/(4*m))*(L-(y0+(h/2)*K_v1)) + (ea/m)*math.atanh(-.9999999999999999) + g
    elif ((u0+(h/2))-(y0+(h/2)*K_v1))/(H-L) >= 1:
        K_v2 = -(k/(4*m))*(L-(y0+(h/2)*K_v1)) + (ea/m)*math.atanh(.9999999999999999) + g
    else:
        K_v2 = -(k/(4*m))*(L-(y0+(h/2)*K_v1)) + (ea/m)*math.atanh(((u0+(h/2))-(y0+(h/2)*K_v1))/(H-L)) + g
   
    if ((u0+(h/2))-(y0+(h/2)*K_v2))/(H-L) <= -1:
        K_v3 = -(k/(4*m))*(L-(y0+(h/2)*K_v2)) + (ea/m)*math.atanh(-.9999999999999999) + g
    elif ((u0+(h/2))-(y0+(h/2)*K_v2))/(H-L) >= 1:
        K_v3 = -(k/(4*m))*(L-(y0+(h/2)*K_v2)) + (ea/m)*math.atanh(.9999999999999999) + g
    else:
        K_v3 = -(k/(4*m))*(L-(y0+(h/2)*K_v2)) + (ea/m)*math.atanh(((u0+(h/2))-(y0+(h/2)*K_v2))/(H-L)) + g
    
    if ((u0+h)-(y0+h*K_v3))/(H-L) <= -1:
        K_v4 = -(k/(4*m))*(L-(y0+h*K_v3)) + (ea/m)*math.atanh(-.9999999999999999) + g
    elif ((u0+h)-(y0+h*K_v3))/(H-L) >= 1:
        K_v4 = -(k/(4*m))*(L-(y0+h*K_v3)) + (ea/m)*math.atanh(.9999999999999999) + g
    else:
        K_v4 = -(k/(4*m))*(L-(y0+h*K_v3)) + (ea/m)*math.atanh(((u0+h)-(y0+h*K_v3))/(H-L)) + g
    
    v1 = v0 + (h/6)*(K_v1 + 2*K_v2 + 2*K_v3 + K_v4)
    
    return y1,v1