"""Katharin Jensen
This code finds the force vs. displacement curve for the compression of a system of parallel units consisting of two
connected rigid links in series with a vertical nonlinear spring, with a linear horizontal spring attached at the 
joint of the two links. 
If desired, it can then create an animation for a small number of links. It can also save this animation as a GIF.
This code can analyze both quasistatic and dynamic motion, with or without damage effects.
For large numbers of links, binning can be applied to reduce computation time"""




#-----------------imports----------------------------------------------------------


import math
from random import random
import time
from numpy import mean
from numpy import linspace
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.animation as animation
import MyRootFinding
import binning
import ani
import dfunc
import RK4Dyn




#-------------initialization-------------------------------------------------------


#begin timing code
start_time = time.clock()


#set type of simulation to run.
dyn = True # sets whether the simulation is dynamic (True) or quasistatic (False).
damage = False #sets whether or not to add damage effects
plasticity = False #sets whether or not to do plasticity with a Jenkins element
dashpot_par = False #sets whether or not to add a dashpot in parallel with the side spring
dashpot_ser = True #sets whether or not to add a dashpot in series with the side spring
bing = False #sets whether or not to use binning
testing = True #sets whether to run in testing mode, with a non-random initial number distribution
auto_dfix = True #if true, if a displacement value is >= H, it will automatically be set to .9999*H. If false, the code won't run if a d value >= H.
multirun = False #sets whether the simulation will run multiple times to make a comparison plot (if true, won't plot average or run animation)


#set animation options
make_ani = True
save_ani = False
auto_frameskip = True #automatically sets frameskip so that if save_ani is on, the animation will play in the same length of time but at 10fps


#set parameters

#basics
ea_t = 10.0 #top spring constant (N)
k_t = 100.0 #side spring constant (N/m)
H = 1.0 #total height (m)
Wm = 4.0 #Weibull modulus
n_l = 5 #number of links
steps = 1000 #number of displacement steps (note, if dyn=True this will be overwritten by t_stop/dt)

#for binning
i_max = 10 #number of bins

#for damage/plasticity
theta_crit = 30.0 #breaking angle, measured from vertical (degrees)
F_slide = 10.0 #critical slipping force for plasticity element (N)

#for the dashpots
c_dash = 1.0

#for dynamics
m = 1.0 #mass (kg)
g = 0.0 #acceleration of gravity (m/s^2)
dt = .02 #time step (s)
t_stop = 10.0 #total time (s)
#if the animation is dynamic, set the number of steps based on the stop time and time step
if dyn:
    steps = int(t_stop/dt)
    print steps

#for animation
xspacing = .5*H #spacing between links in animation (m)
fps = 50 #frame rate (fps) (if save_ani is on, it automatically changes to 10fps, set frameskip accordingly or use auto_frameskip)
frameskip = 1 #frames to skip between displayed frames (e.g. frameskip=3 means play every 3rd frame)
fname = "Example5.gif" #filename for saving GIF (make sure it has the .gif at the end of it)
#automatically set frameskip to play animation at same speed with 10fps if auto_frameskip and save_ani are on
if save_ani and auto_frameskip: 
    #for dynamics
    if dyn:
        frameskip = int((1.0/dt)/10.0)
    else:
        frameskip = int(fps/10)

#for running multiple times
runs = 100 #number of runs for multirun


#create list of displacements using dfunc
#d = dfunc.cospath(steps,3,H,.5)
#d = dfunc.linepath(steps,H)
d1 = dfunc.linepath(steps/2,H*.5)
d2 = dfunc.pausepath(steps/2,H*.5)
#d3 = dfunc.linepath(steps/2,0,d0=H*.5)
#d = dfunc.vpath(steps,1,H,.6)
d=d1+d2



#--------------------check for input errors------------------------------------------------------


#binning related
if bing:
    if i_max >= n_l:
        raise ValueError('The number of bins must be smaller than the number of links')


#damage related
if damage:
    if theta_crit < 0 or theta_crit > 90:
        raise ValueError('The breaking angle must be between 0 and 90 degrees')
    
    
#make sure that only one form of plasticity is on
if (plasticity and (dashpot_par or dashpot_ser)) or (dashpot_par and dashpot_ser):
    raise ValueError('Only ONE of plasticity, dashpot_par, or dashpot_ser may be set to True')



#make sure animation is turned off under conditions that can't be animated

#turn off if there's binning
if bing and make_ani:
    print "Can't create animation if binning is on, turning off animation"
    make_ani = False
    if save_ani:
        save_ani = False

#turn off if there are too many links
if n_l>50 and make_ani:
    print "Too many links for animation, turning off animation"
    make_ani = False
    if save_ani:
        save_ani = False

#turn off if running multiple times
if multirun and make_ani:
    print "Can't create animation when running multiple times, turning off animation"
    make_ani = False
    if save_ani:
        save_ani = False
        
#make sure saving the animation is turned off if the animation isn't being made
if save_ani and not make_ani:
    save_ani = False


#make sure all the numbers that are supposed to be floats are actually floats
#ea_t
if type(ea_t) != float:
    ea_t = float(ea_t)
#k_t    
if type(k_t) != float:
    k_t = float(k_t)
#H    
if type(H) != float:
    H = float(H)
#Wm    
if type(Wm) != float:
    Wm = float(Wm)
#theta_crit    
if type(theta_crit) != float:
    theta_crit = float(theta_crit)
#F_slide    
if type(F_slide) != float:
    F_slide = float(F_slide)
#m    
if type(m) != float:
    m = float(m)
#g    
if type(g) != float:
    g = float(g)
#dt    
if type(dt) != float:
    dt = float(dt)
#t_stop    
if type(t_stop) != float:
    t_stop = float(t_stop)
#xspacing    
if type(xspacing) != float:
    xspacing = float(xspacing)




#---------------------main analysis function-----------------------------------------------------------


def Analysis(ea_t, k_t, H, Wm, n_l, d, i_max, theta_crit, F_slide, c_dash, m, g, dt, t_stop, dyn, damage, plasticity, bing, testing):
    #This function performs the main analysis
    #inputs: all basic, damage, and dynamic parameters, except steps.
    #        list of displacements (d)
    #        all model options
    #outputs: list of force values (P)
    #         list of average force values (P_avg)
    #         list of displacements for the top of the link (y)
    #         critical force for the averate length (F_crit_avg)
    #         list of link lengths (L)
    
    
    #normalize ea and k
    ea = ea_t/n_l
    k = k_t/n_l
    
    
    #generate random numbers
    
    if testing:
        #for testing purposes
        R = [0.1406236293192208, 0.557452455836349, 0.4018884612118805, 0.8610494090625574, 0.005928894753714831]
        #R = linspace(.3*H,.8*H,n_l)
        #R = [0.1406236293192208]
    else:
        #generate list of random numbers
        R = [random() for x in range(n_l)]
    
    
    #create list of link lengths
    
    #function to create length distribution
    def length_distribution(r):
        #apply Weibull distribution
        l = [((math.log(1/num)/math.log(2))**(1.0/Wm)) for num in r]
        #apply exponential distribution to Weibull distribution
        L = [H*(1-math.exp(-num)) for num in l]
        return L
    
    #if binning is applied
    if bing:
        #apply binning to the random number distribution, set number of links to calculate to the new binned number
        W_b,r_b,n = binning.binning(R,i_max)
        #find length distribution
        L = length_distribution(r_b)
        #find average length
        l_avg = length_distribution(R)
        L_avg = mean(l_avg)   
    #if binning is not applied
    else:
        #find length distribution
        L = length_distribution(R)
        #find average length
        L_avg = mean(L)
        #set number of links to calculate to the number of links specified
        n = n_l
    
    
    #initialize lists
    
    #y lists
    y = [[0]*n for x in range(len(d))]
    y_avg = [0]*(len(d))
    
    #P lists
    P = [0]*(len(d))
    P_avg = [0]*(len(d))
    
    #indicators
    #for buckling: 0 = unbuckled, 1 = buckling, 2 = collapsed
    c = [0]*n
    c_avg = 0
    #for breaking: 0 = unbroken, 1 = broken
    b = [0]*n
    b_avg = 0
    
    
    #find critical buckling values
    F_crit = [k*Li/4 for Li in L]
    F_crit_avg = k*L_avg/4
    
    
    #convert theta_crit to radians
    theta_crit = theta_crit*math.pi/180.0
    
    
    #find critical distance for sliding for plasticity
    if plasticity:
        x_crit = F_slide/k
        #initialize distances for plasticity
        x_plas = [0]*n
        x = [0]*n
        x_avg = 0
        x_plas_avg = 0
    elif dashpot_par:
        x = [0]*n
        x_avg = 0
        dx = [0]*n
        dx_avg = 0
    elif dashpot_ser:
        x = [0]*n
        x_avg = 0
        x_dash = [0]*n
        x_dash_avg = 0
        
    
    #set up stuff for root finding if not using dynamic analysis
    if not dyn:
        #if plasticity is not applied 
        if not plasticity:
            #set initial root guesses
            F_old = F_crit #initial guesses
            F_old_avg = F_crit_avg #initial guess for average
            
            #define function to find root of
            def root_g(f,d,L,k,H,ea):
                return L*(1-(4*f)/(k*L)) + (H-L)*math.tanh(f/ea) - d
            
            #define derivative of function
            def root_gprime(f,d,L,k,H,ea):
                return -4/k + ((H-L)/ea)*(1-(math.tanh(f/ea))**2)
            
        #if placticity is applied    
        else:
            #define function to find root of
            def root_g(y,d,L,k,H,ea,x_plas):
                return -ea*math.atanh((d-y)/(H-L)) + k*(L-y)*(.25 - x_plas/(2*math.sqrt(2*L*y-y**2)))
            
            
            
    #set up stuff for ODE solving if using dynamic analysis
    else:
        #initialize velocity
        v = [0]*n
        v_avg = 0
    
    #calculate force at each displacement value
    for i in range(len(d)):
#         if its >= 5:
#             break
        print i   
        
        #initialize list of link forces for this displacement
        F = [0]*n
        
        
        #for each link-------------------------------------------------------------------------------
        for j in range(n):
            #print j    
            
            #if the link isn't buckled, check to see if the displacement has jumped past the top of the link
            if c[j] == 0 and d[i] >= H - L[j]:
                #if so, change the link to buckling
                c[j] = 1 
                
            #if the link isn't buckled
            if c[j] == 0:#-----------------------------------------------------
                #print 'c', 0    
                #calculate the force in the link
                F[j] = ea*math.atanh(d[i]/(H-L[j]))
                #set y
                y[i][j] = 0
                
                #check if the link has buckled
                #if not broken
                if b[j] == 0:
                    #for dynamics
                    if F[j] + m*g >= F_crit[j] and dyn:
                        
                        #set indicator to buckling
                        c[j] = 1
                    #for quasistatic
                    elif F[j] >= F_crit[j] and (not dyn):
                        
                        #set indicator to buckling
                        c[j] = 1
                #if broken
                else:
                    #for dynamics
                    if F[j] + m*g >= 0 and dyn:
                        #print F[j]
                        #set indicator to buckling
                        c[j] = 1
                    #for quasistatic
                    elif F[j] >= 0 and (not dyn):
                        #print F[j]
                        #set indicator to buckling
                        c[j] = 1
                    
            #if the link is buckling
            elif c[j] == 1:#--------------------------------------------------------
                #print 'c', 1 
                #calculate the force in the link
                #for quasistatic
                if not dyn:
                    #if the link isn't broken
                    if b[j] == 0:
                    
                        #set zerotol
                        if L[j] < .3:
                            ztol = .1
                        else:
                            ztol = .01
                        
                            
                        #find the root
                        if not plasticity:
                            #set function arguments
                            args = (d[i],L[j],k,H,ea)
                            #find root
                            F[j] = MyRootFinding.newton_raphson(root_g,F_old[j],root_gprime,args=args,tol=.001,Imax=10*len(d),zerotol=ztol)
                        else:
                            #find spring displacement
                            d_s = x[j] - x_plas[j]
                            #set function arguments
                            args = (d[i],L[j],k,H,ea,x_plas[j])
                            #print 'args', args
                            #find boundaries for bisection method
                            #set lower boundary to either 0 or the smallest number for atanh
                            if d[i]-H+L[j] < 0:
                                b_low = .0001
                            else:
                                b_low = (d[i]-H+L[j])*1.0001
                            #set upper boundary to the largest number for atanh
                            b_high = (H-L[j]+d[i])*.9999
                            y[i][j] = MyRootFinding.bisection(root_g,b_low,b_high,args=args)
                        
                        #Find other required values
                        if plasticity:
                            #store y value
                            #find force
                            F[j] = ea*math.atanh((d[i]-y[i][j])/(H-L[j]))
                            #find new x
                            x[j] = math.sqrt(2*L[j]*y[i][j]-y[i][j]**2)/2
                            #find new x_plas if x_plas is changing
                            if d_s > x_crit:
                                x_plas[j] = x[j] - x_crit
                            elif d_s < -x_crit:
                                x_plas[j] = x[j] + x_crit
                            #print 'x_plas', x_plas[j]
                        else:
                            #assign current value as next guess value, to be more accurate and help prevent diverging
                            F_old[j] = F[j]
                            #find and store y value
                            y[i][j] = L[j] - 4*F[j]/k
                        
                    #if the link is broken
                    else:
                        #the side spring force is 0, so the applied force is also 0
                        F[j] = 0
                        #set y value
                        y[i][j] = d[i]
                    
                #for dynamic
                else:
                    #solve ODE for y with 4th order Runge Kutta method
                    
                    #if the link isn't broken
                    if b[j] == 0:
                        if dashpot_par:
                            #save previous x value
                            xold = x[j]
                            #find final argument value
                            fin_arg = -c_dash*dx[j]/k
                            #run Runge Kutta method for plasticity with fin_arg instead of x_plas to find y and the velocity
                            y[i][j],v[j] = RK4Dyn.RK4_plaslink(dt,d[i-1],d[i],y[i-1][j],v[j],H,L[j],k,ea,m,g,fin_arg)
                            #find new x
                            if y[i][j] < 0:
                                x[j] = 0
                            elif y[i][j] > L[j]:
                                x[j] = L[j]/2
                            else:
                                x[j] = math.sqrt(2*L[j]*y[i][j]-y[i][j]**2)/2
                            #find new dx
                            dx[j] = (x[j] - xold)/dt
                            
                        elif dashpot_ser:
                            #run Runge Kutta method for plasticity with x_dash instead of x_plas to find y and the velocity
                            y[i][j],v[j] = RK4Dyn.RK4_plaslink(dt,d[i-1],d[i],y[i-1][j],v[j],H,L[j],k,ea,m,g,x_dash[j])
                            #get new x
                            if y[i][j] < 0:
                                x[j] = 0
                            elif y[i][j] > L[j]:
                                x[j] = L[j]/2
                            else:
                                x[j] = math.sqrt(2*L[j]*y[i][j]-y[i][j]**2)/2
                            #use Runge Kutta method to get new x_dash
                            x_dash[j] = RK4Dyn.RK4_dashser(x_dash[j],dt,k,x[j],c_dash)
                            
                        elif plasticity:
                            #find spring displacement
                            d_s = x[j] - x_plas[j]
                            #run Runge Kutta method for plasticity to find y and the velocity
                            y[i][j],v[j] = RK4Dyn.RK4_plaslink(dt,d[i-1],d[i],y[i-1][j],v[j],H,L[j],k,ea,m,g,x_plas[j])
                            #find new x
                            if y[i][j] < 0:
                                x[j] = 0
                            elif y[i][j] > L[j]:
                                x[j] = L[j]/2
                            else:
                                x[j] = math.sqrt(2*L[j]*y[i][j]-y[i][j]**2)/2
                            #find new x_plas if x_plas is changing
                            if d_s > x_crit:
                                x_plas[j] = x[j] - x_crit
                            elif d_s < -x_crit:
                                x_plas[j] = x[j] + x_crit
                        else:
                            #run regular Runge Kutta method to find y and the velocity
                            y[i][j],v[j] = RK4Dyn.RK4_link(dt,d[i-1],d[i],y[i-1][j],v[j],H,L[j],k,ea,m,g)
                        
                        
                    else:
                        #run Runge Kutta method to find y and the velocity, with k=0 for breakage
                        y[i][j],v[j] = RK4Dyn.RK4_link(dt,d[i-1],d[i],y[i-1][j],v[j],H,L[j],0,ea,m,g)
                    #correct v value if necessary to prevent divergence
                    if v[j] > 1000 or math.isnan(v[j]):
                        v[j] = 1000
                    
                    #correct y value if necessary
                    #make sure the displacement doesn't pass the top of the link
                    if d[i]-y[i][j] >= H-L[j]:
                        y[i][j] = d[i] - (H-L[j])*.999999
                    #make sure the spring doesn't stretch too long
                    if d[i]-y[i][j] <= L[j]-H:
                        y[i][j] = d[i] - (L[j]-H)*.999999
                    
                    #calculate force
                    F[j] = ea*math.atanh((d[i]-y[i][j])/(H-L[j]))
                
                #if damage is turned on and the link hasn't broken, check to see if it has
                if b[j] == 0 and damage:
                    #calculate theta
                    #if y steps below zero, making the argument for acos > 1, make sure theta is zero
                    if y[i][j] < 0:
                        theta = 0
                    else:
                        theta = math.acos((L[j]-y[i][j])/L[j])
                    #check if theta is larger than the breaking angle
                    if theta >= theta_crit:
                        #set to broken
                        b[j] = 1
                
                #check to see if the link has completely collapsed
                if y[i][j] >= L[j]:
                    #set y correctly
                    y[i][j] = L[j]
                    
                    #set indicator to collapsed
                    c[j] = 2
                    #recalculate force to reflect collapse
                    F[j] = ea*math.atanh((d[i] - L[j])/(H - L[j]))
                    #set velocity to zero if running dynamics
                    if dyn:
                        v[j] = 0
                
                #check to see if the link has unbuckled
                if not plasticity:
                    if y[i][j] <= 0:
                        #set y correctly
                        y[i][j] = 0
                        #set indicator to unbuckled
                        c[j] = 0
                        #recalculate force to reflect unbuckling
                        #if the atanh argument is >= 1
                        if d[i] >= (H-L[j]):
                            #make the force as large as it can be
                            F[j] = ea*math.atanh(.9999999999999999)
                        #otherwise calculate the force normally
                        else:
                            F[j] = ea*math.atanh(d[i]/(H - L[j]))
                else:
                    if (y[i][j] > 10*y[i-1][j] and y[i-1][j] != 0 and not dyn) or (y[i][j] < 0 and dyn):
                        #set y correctly
                        y[i][j] = 0
                        #set indicator to unbuckled
                        c[j] = 0
                        #recalculate force to reflect unbuckling
                        #if the atanh argument is >= 1
                        if d[i] >= (H-L[j]):
                            #make the force as large as it can be
                            F[j] = ea*math.atanh(.9999999999999999)
                        #otherwise calculate the force normally
                        else:
                            F[j] = ea*math.atanh(d[i]/(H - L[j]))
            
            #if the link is completely collapsed
            elif c[j] == 2:#------------------------------------------------------------------
                #print 'c', 2 
                #store y
                y[i][j] = L[j]

                #calculate the force in the link
                #if the atanh argument <= -1
                if (d[i] - L[j])/(H - L[j]) <= -1:
                    #make the force as small as possible
                    F[j] = ea*math.atanh(-.9999999999999999)
                #otherwise calculate normally
                else:
                    F[j] = ea*math.atanh((d[i] - L[j])/(H - L[j]))
                
                #check to see if the link is uncollapsing
                #for dynamics, check if the spring force will overcome gravity
                if F[j] < -m*g and dyn:
                    #set indicator to buckling
                    c[j] = 1
                #for quasistatic, check if the spring force < 0
                elif F[j] < 0 and (not dyn):
                    #set indicator to buckling
                    c[j] = 1
                    
        
        #calculate total force
        
        #if binning is applied
        if bing:
            P[i] = 0 #initialize force sum
            for j in range(n):
                P[i] += F[j]*W_b[j]*n_l #sum weighted forces
        #if binning is not applied        
        else:    
            P[i] = sum(F)
            
            
            
        #for the average force------------------------------------------------------------------------------------
        
        
        #if the link isn't buckled, check to see if the displacement has jumped past the top of the link
        if c_avg == 0 and d[i] >= H - L_avg:
            #if so, change the link to buckling
            c_avg = 1
            
            
        #if the link isn't buckled
        if c_avg == 0:#-----------------------------------------------------
            
            #calculate the force in the link
            F_avg = ea*math.atanh(d[i]/(H-L_avg))
            #set y
            y_avg[i] = 0
            
            #check if the link has buckled
            #for dynamics
            if F_avg + m*g >= F_crit_avg and dyn:
                #set indicator to buckling
                c_avg = 1
            #for quasistatic
            elif F_avg >= F_crit_avg and (not dyn):
                #set indicator to buckling
                c_avg = 1
                
                
        #if the link is buckling
        elif c_avg == 1:#---------------------------------------------------
            
            #calculate the force in the link
            #for quasistatic
            if not dyn:
                #if the link isn't broken
                if b_avg == 0:
                
                    #set the guess
                    if plasticity:
                        if y_avg[i-1] == 0:
                            guess = .00001
                        else:
                            guess = y_avg[i-1]
                    else:
                        guess = F_old_avg
                    
                    #find the root
                    if not plasticity:
                            #set function arguments
                            args = (d[i],L_avg,k,H,ea)
                            #find root
                            F_avg = MyRootFinding.newton_raphson(root_g,guess,root_gprime,args=(d[i],L_avg,k,H,ea),tol=.001,Imax=10*len(d),zerotol=.01)
                    else:
                        #find spring displacement
                        d_s_avg = x_avg - x_plas_avg
                        #set function arguments
                        args_avg = (d[i],L_avg,k,H,ea,x_plas_avg)
                        #print 'args', args
                        #find boundaries for bisection method
                        #set lower boundary to either 0 or the smallest number for atanh
                        if d[i]-H+L_avg < 0:
                            b_low_avg = .0001
                        else:
                            b_low_avg = (d[i]-H+L_avg)*1.0001
                        #set upper boundary to the largest number for atanh
                        b_high_avg = (H-L_avg+d[i])*.9999
                        y_avg[i] = MyRootFinding.bisection(root_g,b_low_avg,b_high_avg,args=args_avg)
                        
                    #Find other required values
                    if plasticity:
                        #store y value
                        #find force
                        F_avg = ea*math.atanh((d[i]-y_avg[i])/(H-L_avg))
                        print 'F', F_avg
                        #find new x
                        x_avg = math.sqrt(2*L_avg*y_avg[i]-y_avg[i]**2)/2
                        #find new x_plas if x_plas is changing
                        if d_s_avg > x_crit:
                            x_plas_avg = x_avg - x_crit
                        elif d_s_avg < -x_crit:
                            x_plas_avg = x_avg + x_crit
                        #print 'x_plas', x_plas[j]
                    else:
                        #assign current value as next guess value, to be more accurate and help prevent diverging
                        F_old_avg = F_avg
                        #find and store y value
                        y_avg[i] = L_avg - 4*F_avg/k
                        
                        
                #if the link is broken
                else:
                    #the side spring force is 0, so the applied force is also 0
                    F_avg = 0
                    #set y value
                    y_avg[i] = d[i]
                
            #for dynamic
            else:
                #solve ODE for y with 4th order Runge Kutta method
                    
                #if the link isn't broken
                if b[j] == 0:
                    if dashpot_par:
                        #save previous x value
                        xold_avg = x_avg
                        #find final argument value
                        fin_arg_avg = -c_dash*dx_avg/k
                        #run Runge Kutta method for plasticity with fin_arg instead of x_plas to find y and the velocity
                        y_avg[i],v_avg = RK4Dyn.RK4_plaslink(dt,d[i-1],d[i],y_avg[i-1],v_avg,H,L_avg,k,ea,m,g,fin_arg_avg)
                        #find new x
                        if y_avg[i] < 0:
                            x_avg = 0
                        elif y_avg[i] > L_avg:
                            x_avg = L_avg/2
                        else:
                            x_avg = math.sqrt(2*L_avg*y_avg[i]-y_avg[i]**2)/2
                        #find new dx
                        dx_avg = (x_avg - xold_avg)/dt
                        
                    elif dashpot_ser:
                        #run Runge Kutta method for plasticity with x_dash instead of x_plas to find y and the velocity
                        y_avg[i],v_avg = RK4Dyn.RK4_plaslink(dt,d[i-1],d[i],y_avg[i-1],v_avg,H,L_avg,k,ea,m,g,x_dash_avg)
                        #get new x
                        if y_avg[i] < 0:
                            x_avg = 0
                        elif y_avg[i] > L_avg:
                            x_avg = L_avg/2
                        else:
                            x_avg = math.sqrt(2*L_avg*y_avg[i]-y_avg[i]**2)/2
                        #use Runge Kutta method to get new x_dash
                        x_dash_avg = RK4Dyn.RK4_dashser(x_dash_avg,dt,k,x_avg,c_dash)
                        
                    elif plasticity:
                        #find spring displacement
                        d_s_avg = x_avg - x_plas_avg
                        #run Runge Kutta method for plasticity to find y and the velocity
                        y_avg[i],v_avg = RK4Dyn.RK4_plaslink(dt,d[i-1],d[i],y_avg[i-1],v_avg,H,L_avg,k,ea,m,g,x_plas_avg)
                        #find new x
                        if y_avg[i] < 0:
                            x_avg = 0
                        elif y_avg[i] > L_avg:
                            x_avg = L_avg/2
                        else:
                            x_avg = math.sqrt(2*L_avg*y_avg[i]-y_avg[i]**2)/2
                        #find new x_plas if x_plas is changing
                        if d_s_avg > x_crit:
                            x_plas_avg = x_avg - x_crit
                        elif d_s_avg < -x_crit:
                            x_plas_avg = x_avg + x_crit
                    else:
                        #run Runge Kutta method to find y and the velocity
                        y_avg[i],v_avg = RK4Dyn.RK4_link(dt,d[i-1],d[i],y_avg[i-1],v_avg,H,L_avg,k,ea,m,g)

                else:
                    #run Runge Kutta method to find y and the velocity, with k=0 for breakage
                    y_avg[i],v_avg = RK4Dyn.RK4_link(dt,d[i-1],d[i],y_avg[i-1],v_avg,H,L_avg,0,ea,m,g)

                #correct v value if necessary to prevent divergence
                if v_avg > 1000 or math.isnan(v_avg):
                    v_avg = 1000
                
                #correct y value if necessary
                #make sure the displacement doesn't pass the top of the link
                if d[i]-y_avg[i] >= H-L_avg:
                    y_avg[i] = d[i] - (H-L_avg)*.999999
                #make sure the spring doesn't stretch too long
                if d[i]-y_avg[i] <= L_avg-H:
                    y_avg[i] = d[i] - (L_avg-H)*.999999
            
                #calculate force using the calculated y value
                F_avg = ea*math.atanh((d[i]-y_avg[i])/(H-L_avg))
            
            #if damage is on and the link hasn't broken, check to see if it has
            if b_avg == 0 and damage:
                #calculate theta
                #if y steps below zero, making the argument for acos > 1, make sure theta is zero
                if y_avg[i] < 0:
                    theta = 0
                else:
                    theta = math.acos((L_avg-y_avg[i])/L_avg)
                #check if theta is larger than the breaking angle
                if theta >= theta_crit:
                    #set to broken
                    b_avg = 1
            
            #check to see if the link has completely collapsed
            if y_avg[i] >= L_avg:
                #set y correctly
                y_avg[i] = L_avg
                #set indicator to collapsed
                c_avg = 2
                #recalculate force to reflect collapse
                F_avg = ea*math.atanh((d[i] - L_avg)/(H - L_avg))
                #set velocity to zero if running dynamics
                if dyn:
                    v_avg = 0
            
            #check to see if the link has unbuckled
            if not plasticity:
                if y_avg[i] < 0:
                    #set y correctly
                    y_avg[i] = 0
                    #set indicator to unbuckled
                    c_avg = 0
                    #recalculate force to reflect unbuckling
                    F_avg = ea*math.atanh(d[i]/(H - L_avg))
            else:
                if (y_avg[i] > 10*y_avg[i-1] and y_avg[i-1] != 0 and not dyn) or (y_avg[i] < 0 and dyn):
                    #set y correctly
                    y_avg[i] = 0
                    #set indicator to unbuckled
                    c_avg = 0
                    #recalculate force to reflect unbuckling
                    #if the atanh argument is >= 1
                    if d[i] >= (H-L_avg):
                        #make the force as large as it can be
                        F_avg = ea*math.atanh(.9999999999999999)
                    #otherwise calculate the force normally
                    else:
                        F_avg = ea*math.atanh(d[i]/(H - L_avg))
        
        
        #if the link is completely collapsed
        elif c_avg == 2:#-----------------------------------------------------------
            
            #store y
            y_avg[i] = L_avg
            
            #calculate the force in the link
            #if the atanh argument is <= -1
            if (d[i] - L_avg)/(H - L_avg) <= -1:
                #make force as small as possible
                F_avg = ea*math.atanh(-.9999999999999999)
            #otherwise calculate force normally
            else:
                F_avg = ea*math.atanh((d[i] - L_avg)/(H - L_avg))
            
            #check to see if the link is uncollapsing
            #for dynamics, check if the spring force will overcome gravity
            if F_avg < -m*g and dyn:
                #set indicator to buckling
                c_avg = 1
            #for quasistatic, check if the spring force < 0
            elif F_avg < 0 and (not dyn):
                #set indicator to buckling
                c_avg = 1
                    
        
        #calculate total average force
        P_avg[i] = F_avg*n_l
    
    
    
    
    return P, P_avg, y, F_crit_avg, L




#-------------------------------run simulation----------------------------------------------------------------------


#make sure none of the displacements are >= H
#check each displacement value
for i in range(len(d)):
    #if the displacement >= H
    if d[i] >= H:
        #if auto_dfix is on, set d[i] to .9999*H
        if auto_dfix:
            d[i] = .9999*H
        #if auto displacement fixing is off
        else:
            #raise error
            raise ValueError("Found d >= H. Change d so that all d < H and try again.")


#run analysis function
 
#if multirun is on
if multirun:
    #initialize force lists
    P = [[] for i in range(runs)]
    F_crit_avg = [0]*runs
    #run the simulation multiple times
    for k in range(runs):
        print "run", k
        #find force and F_crit_avg for each run
        P[k], P_avg, y, F_crit_avg[k], L = Analysis(ea_t, k_t, H, Wm, n_l, d, i_max, theta_crit, F_slide, c_dash, m, g, dt, t_stop, dyn, damage, plasticity, bing, testing)
#if only running once
else:
    P, P_avg, y, F_crit_avg, L = Analysis(ea_t, k_t, H, Wm, n_l, d, i_max, theta_crit, F_slide, c_dash, m, g, dt, t_stop, dyn, damage, plasticity, bing, testing)
    

#stop timing code
print "Elapsed time is %f seconds" % (time.clock() - start_time)

    
#plot

#if multirun is on
if multirun:
    for k in range(runs):
        plt.plot(d,P[k],'b',lw=2,alpha=0.2)
    #set axes
    if not dyn:    
        plt.axis([0,H+.1*H,0,max(F_crit_avg)*n_l*1.5])
    else:
        #find min and max force values
        #initialize
        minP = 0
        maxP = 0
        for k in range(runs):
            minPk = min(P[k])
            maxPk = max(P[k])
            if minPk < minP:
                minP = minPk
            if maxPk > maxP:
                maxP = maxPk
        plt.axis([0,H+.1*H,minP*1.1,maxP*1.1])
    plt.xlabel("Displacement (m)")
    plt.ylabel("Force (N)")
    plt.title("%d %d-link Force vs. Displacement Curves" % (runs,n_l))
#if running once
else:
    plt.plot(d,P,d,P_avg,'r--',lw=1)
    #set the axis
    #for quasistatic:
    if not dyn:
        plt.axis([0,H+.1*H,0,F_crit_avg*n_l*1.1])
    #for dynamics
    else:
        #plt.axis([0,H+.1*H,0,20])
        plt.axis([0,H+.1*H,min(P_avg)*1.1,max(P_avg)*1.1])
    plt.xlabel("Displacement (m)")
    plt.ylabel("Force (N)")
    plt.title("Force vs. Displacement")
    plt.legend(["Actual Force Value","Force Value for Average Length"],fontsize = 'medium')



#show plot
plt.show()




#--------------------------animate----------------------------------------------------------------------------


if make_ani:
    #create animation
    anim = ani.ani(L,y,H,n_l,d,P,P_avg,xspacing,F_crit_avg,save_ani,fps=fps,frameskip=frameskip,dyn=dyn,t_step=dt)

if save_ani:
    #save animation as GIF
    start_ani = time.clock()
    print "begin saving animation"
    anim.save(fname, writer='imagemagick')
    print "done saving animation"
    print "animation save time is %f seconds" % (time.clock() - start_ani)

