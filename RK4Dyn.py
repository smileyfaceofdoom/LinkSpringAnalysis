"""Katharin Jensen
This code implements the 4th order Runge-Kutta method for the dynamic analysis equations"""

#imports
import math


def RK4(dt,u0,u1,y0,v0,H,L,k,ea,m,g):
    #inputs: u at the beginning (u0) and end (u1) of the step
    #        y and v at the beginning of the time step
    #        parameters H, L, k, ea, and m
    #outputs: y and v at the end of the time step
    
    
    #find displacement step size du
    du = u1-u0
    
    
    #find time step size h
    h = dt
    
    
    #find next v and y values
    
    #Find K_v1
    if ((u0)-(y0))/(H-L) <= -1:
        K_v1 = (k/(4*m))*(L-(y0)) - (ea/m)*math.atanh(-.9999999999999999) - g
        
    elif ((u0)-(y0))/(H-L) >= 1:
        K_v1 = (k/(4*m))*(L-(y0)) - (ea/m)*math.atanh(.9999999999999999) - g
        
    else:
        K_v1 = (k/(4*m))*(L-(y0)) - (ea/m)*math.atanh(((u0)-(y0))/(H-L)) - g
    
    #find K_y1
    K_y1 = -v0
    
    
    #find K_v2
    if ((u0+(du/2))-(y0+(h/2)*K_y1))/(H-L) <= -1:
        K_v2 = (k/(4*m))*(L-(y0+(h/2)*K_y1)) - (ea/m)*math.atanh(-.9999999999999999) - g
        
    elif ((u0+(du/2))-(y0+(h/2)*K_y1))/(H-L) >= 1:
        K_v2 = (k/(4*m))*(L-(y0+(h/2)*K_y1)) - (ea/m)*math.atanh(.9999999999999999) - g
        
    else:
        K_v2 = (k/(4*m))*(L-(y0+(h/2)*K_y1)) - (ea/m)*math.atanh(((u0+(du/2))-(y0+(h/2)*K_y1))/(H-L)) - g
    
    #find K_y2
    K_y2 = -(v0 + (h/2)*K_v1)
    
    
    #find K_v3
    if ((u0+(du/2))-(y0+(h/2)*K_y2))/(H-L) <= -1:
        K_v3 = (k/(4*m))*(L-(y0+(h/2)*K_y2)) - (ea/m)*math.atanh(-.9999999999999999) - g
        
    elif ((u0+(du/2))-(y0+(h/2)*K_y2))/(H-L) >= 1:
        K_v3 = (k/(4*m))*(L-(y0+(h/2)*K_y2)) - (ea/m)*math.atanh(.9999999999999999) - g
        
    else:
        K_v3 = (k/(4*m))*(L-(y0+(h/2)*K_y2)) - (ea/m)*math.atanh(((u0+(du/2))-(y0+(h/2)*K_y2))/(H-L)) - g
        
    #find K_y3
    K_y3 = -(v0 + (h/2)*K_v2)
    
    
    #find K_v4
    if ((u1)-(y0+h*K_y3))/(H-L) <= -1:
        K_v4 = (k/(4*m))*(L-(y0+h*K_y3)) - (ea/m)*math.atanh(-.9999999999999999) - g
        
    elif ((u1)-(y0+h*K_y3))/(H-L) >= 1:
        K_v4 = (k/(4*m))*(L-(y0+h*K_y3)) - (ea/m)*math.atanh(.9999999999999999) - g
        
    else:
        K_v4 = (k/(4*m))*(L-(y0+h*K_y3)) - (ea/m)*math.atanh(((u1)-(y0+h*K_y3))/(H-L)) - g
    
    #find K_y4
    K_y4 = -(v0 + h*K_v3)
        
    
    #find next v value
    v1 = v0 + (h/6)*(K_v1 + 2*K_v2 + 2*K_v3 + K_v4)
    
    #find next y value
    y1 = y0 + (h/6)*(K_y1 + 2*K_y2 + 2*K_y3 + K_y4)
    
    return y1,v1