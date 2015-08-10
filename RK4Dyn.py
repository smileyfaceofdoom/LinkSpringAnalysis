"""Katharin Jensen
This code implements the 4th order Runge-Kutta method for the dynamic analysis equations"""

#imports
import math


def RK4_point(dt,u0,u1,y0,v0,H,L,k,ea,m,g):
    #This function implements the 4th order Runge Kutta method for the point mass dynamics
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



def RK4_link(dt,u0,u1,y0,v0,H,L,k,ea,m,g):
    #This function implements the 4th order Runge Kutta method for the link mass dynamics
    #inputs: dt, the size of the time step
    #        u at the beginning (u0) and end (u1) of the step
    #        y and v at the beginning of the time step
    #        parameters H, L, k, ea, m, and g
    #outputs: y and v at the end of the time step
    
    
    #find displacement step size du
    du = u1 - u0
    
    if y0 == 0:
        y0 = .0001
#         v0 = 0
    
    #get time step size h
    h = dt

    
    #find K_v1
    #prevent atanh argument from going out of bounds
    if ((u0)-(y0))/(H-L) <= -1:
        iht1 = -.9999999999999999
    elif ((u0)-(y0))/(H-L) >= 1:
        iht1 = .9999999999999999
    else:
        iht1 = ((u0)-(y0))/(H-L)
    
    K_v1 = (-((v0)**2*(L-(y0)))/(2*L*(y0)-(y0)**2) - 
            (2*L*(y0)-(y0)**2)*((m*g/2  + ea*math.atanh(iht1) - k*(L-(y0))/4)/((23./48.)*m*L**2 + m*L*(y0) - .5*m*(y0)**2)))
    
    #find K_y1
    K_y1 = -(v0)
    
    #find K_v2
    #prevent atanh argument from going out of bounds
    if ((u0+du/2)-(y0+(h/2)*K_y1))/(H-L) <= -1:
        iht2 = -.9999999999999999
    elif ((u0+du/2)-(y0+(h/2)*K_y1))/(H-L) >= 1:
        iht2 = .9999999999999999
    else:
        iht2 = ((u0+du/2)-(y0+(h/2)*K_y1))/(H-L)
    
    K_v2 = (-((v0+(h/2)*K_v1)**2*(L-(y0+(h/2)*K_y1)))/(2*L*(y0+(h/2)*K_y1)-(y0+(h/2)*K_y1)**2) - 
            (2*L*(y0+(h/2)*K_y1)-(y0+(h/2)*K_y1)**2)*((m*g/2  + ea*math.atanh(iht2) - k*(L-(y0+(h/2)*K_y1))/4)/((23./48.)*m*L**2 + m*L*(y0+(h/2)*K_y1) - .5*m*(y0+(h/2)*K_y1)**2)))
    
    #find K_y2
    K_y2 = -(v0+(h/2)*K_v1)
    
    #find K_v3
    #prevent atanh argument from going out of bounds
    if ((u0+du/2)-(y0+(h/2)*K_y2))/(H-L) <= -1:
        iht3 = -.9999999999999999
    elif ((u0+du/2)-(y0+(h/2)*K_y2))/(H-L) >= 1:
        iht3 = .9999999999999999
    else:
        iht3 = ((u0+du/2)-(y0+(h/2)*K_y2))/(H-L)
    
    K_v3 = (-((v0+(h/2)*K_v2)**2*(L-(y0+(h/2)*K_y2)))/(2*L*(y0+(h/2)*K_y2)-(y0+(h/2)*K_y2)**2) - 
            (2*L*(y0+(h/2)*K_y2)-(y0+(h/2)*K_y2)**2)*((m*g/2  + ea*math.atanh(iht3) - k*(L-(y0+(h/2)*K_y2))/4)/((23./48.)*m*L**2 + m*L*(y0+(h/2)*K_y2) - .5*m*(y0+(h/2)*K_y2)**2)))
    
    #find K_y3
    K_y3 = -(v0+(h/2)*K_v2)
    
    #find K_v4
    #prevent atanh argument from going out of bounds
    if ((u1)-(y0+h*K_y3))/(H-L) <= -1:
        iht4 = -.9999999999999999
    elif ((u1)-(y0+h*K_y3))/(H-L) >= 1:
        iht4 = .9999999999999999
    else:
        iht4 = ((u1)-(y0+h*K_y3))/(H-L)
    
    K_v4 = (-((v0+h*K_v3)**2*(L-(y0+h*K_y3)))/(2*L*(y0+h*K_y3)-(y0+h*K_y3)**2) - 
            (2*L*(y0+h*K_y3)-(y0+h*K_y3)**2)*((m*g/2  + ea*math.atanh(iht4) - k*(L-(y0+h*K_y3))/4)/((23./48.)*m*L**2 + m*L*(y0+h*K_y3) - .5*m*(y0+h*K_y3)**2)))
    
    #find K_y4
    K_y4 = -(v0+h*K_v3)
    
    #find new v
    v1 = v0 + (h/6)*(K_v1 + 2*K_v2 + 2*K_v3 + K_v4)

    
    #find new y
    y1 = y0 + (h/6)*(K_y1 + 2*K_y2 + 2*K_y3 + K_y4)
    
    
    return y1, v1
    



def RK4_plaslink(dt,u0,u1,y0,v0,H,L,k,ea,m,g,x_plas):
    #This function implements the 4th order Runge Kutta method for the link mass dynamics with plasticity
    #inputs: dt, the size of the time step
    #        u at the beginning (u0) and end (u1) of the step
    #        y and v at the beginning of the time step
    #        parameters H, L, k, ea, m, and g
    #outputs: y and v at the end of the time step
    
    
    #find displacement step size du
    du = u1 - u0
    
    if y0 == 0:
        y0 = .0001
#         v0 = 0
    
    #get time step size h
    h = dt

    
    #find K_v1
    #prevent atanh argument from going out of bounds
    if ((u0)-(y0))/(H-L) <= -1:
        iht1 = -.9999999999999999
    elif ((u0)-(y0))/(H-L) >= 1:
        iht1 = .9999999999999999
    else:
        iht1 = ((u0)-(y0))/(H-L)
        
    #prevent square root argument from going out of bounds
    if 2*L*(y0)-(y0)**2 < 0:
        sqt1 = 0
    else:
        sqt1 = 2*L*(y0)-(y0)**2
    
    K_v1 = (-((v0)**2*(L-(y0)))/(2*L*(y0)-(y0)**2) - 
            ((2*L*(y0)-(y0)**2)*(m*g/2.  + ea*math.atanh(iht1) - k*(L-(y0))/4) + k*x_plas*math.sqrt(sqt1)*((L-(y0))/2))
            /((23./48.)*m*L**2 + m*L*(y0) - .5*m*(y0)**2))
    
    #find K_y1
    K_y1 = -(v0)
    
    #find K_v2
    #prevent atanh argument from going out of bounds
    if ((u0+du/2)-(y0+(h/2)*K_y1))/(H-L) <= -1:
        iht2 = -.9999999999999999
    elif ((u0+du/2)-(y0+(h/2)*K_y1))/(H-L) >= 1:
        iht2 = .9999999999999999
    else:
        iht2 = ((u0+du/2)-(y0+(h/2)*K_y1))/(H-L)
    
    #prevent square root argument from going out of bounds
    if 2*L*(y0+(h/2)*K_y1)-(y0+(h/2)*K_y1)**2 < 0:
        sqt2 = 0
    else:
        sqt2 = 2*L*(y0+(h/2)*K_y1)-(y0+(h/2)*K_y1)**2
    
    K_v2 = (-((v0+(h/2)*K_v1)**2*(L-(y0+(h/2)*K_y1)))/(2*L*(y0+(h/2)*K_y1)-(y0+(h/2)*K_y1)**2) - 
            ((2*L*(y0+(h/2)*K_y1)-(y0+(h/2)*K_y1)**2)*(m*g/2  + ea*math.atanh(iht2) - k*(L-(y0+(h/2)*K_y1))/4) + k*x_plas*math.sqrt(sqt2)*((L-(y0+(h/2)*K_y1))/2))
            /((23./48.)*m*L**2 + m*L*(y0+(h/2)*K_y1) - .5*m*(y0+(h/2)*K_y1)**2))
    
    #find K_y2
    K_y2 = -(v0+(h/2)*K_v1)
    
    #find K_v3
    #prevent atanh argument from going out of bounds
    if ((u0+du/2)-(y0+(h/2)*K_y2))/(H-L) <= -1:
        iht3 = -.9999999999999999
    elif ((u0+du/2)-(y0+(h/2)*K_y2))/(H-L) >= 1:
        iht3 = .9999999999999999
    else:
        iht3 = ((u0+du/2)-(y0+(h/2)*K_y2))/(H-L)
    
    #prevent square root argument from going out of bounds
    if 2*L*(y0+(h/2)*K_y2)-(y0+(h/2)*K_y2)**2 < 0:
        sqt3 = 0
    else:
        sqt3 = 2*L*(y0+(h/2)*K_y2)-(y0+(h/2)*K_y2)**2
    
    K_v3 = (-((v0+(h/2)*K_v2)**2*(L-(y0+(h/2)*K_y2)))/(2*L*(y0+(h/2)*K_y2)-(y0+(h/2)*K_y2)**2) - 
            ((2*L*(y0+(h/2)*K_y2)-(y0+(h/2)*K_y2)**2)*(m*g/2  + ea*math.atanh(iht3) - k*(L-(y0+(h/2)*K_y2))/4) + k*x_plas*math.sqrt(sqt3)*((L-(y0+(h/2)*K_y2))/2))
            /((23./48.)*m*L**2 + m*L*(y0+(h/2)*K_y2) - .5*m*(y0+(h/2)*K_y2)**2))
    
    #find K_y3
    K_y3 = -(v0+(h/2)*K_v2)
    
    #find K_v4
    #prevent atanh argument from going out of bounds
    if ((u1)-(y0+h*K_y3))/(H-L) <= -1:
        iht4 = -.9999999999999999
    elif ((u1)-(y0+h*K_y3))/(H-L) >= 1:
        iht4 = .9999999999999999
    else:
        iht4 = ((u1)-(y0+h*K_y3))/(H-L)
    
    #prevent square root argument from going out of bounds
    if 2*L*(y0+h*K_y3)-(y0+h*K_y3)**2 < 0:
        sqt4 = 0
    else:
        sqt4 = 2*L*(y0+h*K_y3)-(y0+h*K_y3)**2
    
    K_v4 = (-((v0+h*K_v3)**2*(L-(y0+h*K_y3)))/(2*L*(y0+h*K_y3)-(y0+h*K_y3)**2) - 
            ((2*L*(y0+h*K_y3)-(y0+h*K_y3)**2)*(m*g/2  + ea*math.atanh(iht4) - k*(L-(y0+h*K_y3))/4) + k*x_plas*math.sqrt(sqt4)*((L-(y0+h*K_y3))/2))
            /((23./48.)*m*L**2 + m*L*(y0+h*K_y3) - .5*m*(y0+h*K_y3)**2))
    
    #find K_y4
    K_y4 = -(v0+h*K_v3)
    
    #find new v
    v1 = v0 + (h/6)*(K_v1 + 2*K_v2 + 2*K_v3 + K_v4)

    
    #find new y
    y1 = y0 + (h/6)*(K_y1 + 2*K_y2 + 2*K_y3 + K_y4)
    
    
    return y1, v1




