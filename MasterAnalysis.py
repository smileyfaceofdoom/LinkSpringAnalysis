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




#-------------initialization-------------------------------------------------------


#begin timing code
start_time = time.clock()


#set type of simulation to run.
dyn = True # sets whether the simulation is dynamic (True) or quasistatic (False).
damage = True #sets whether or not to add damage effects
bing = False #sets whether or not to use binning
testing = True #sets whether to run in testing mode, with a non-random initial number distribution
auto_dfix = True #If true, if a displacement value is >= H, it will automatically be set to .9999*H. If false, the code won't run if a d value >= H.
multirun = False #sets whether the simulation will run multiple times to make a comparison plot (if true, won't plot average or run animation)


#set animation options
make_ani = True
save_ani = False
auto_frameskip = True #automatically sets frameskip so that if save_ani is on, the animation will play in the same length of time but at 10fps


#set parameters

#basics
ea_t = 10.0 #top spring constant (N)
k_t = 100.0 #side spring constant (N/m)
H = 2.0 #total height (m)
Wm = 4.0 #Weibull modulus
n_l = 4 #number of links
steps = 1000 #number of displacement steps (note, if dyn=True this will be overwritten by t_stop/dt)

#for binning
i_max = 21 #number of bins

#for damage
theta_crit = 30.0 #breaking angle, measured from vertical (degrees)

#for dynamics
m = 1.0 #mass (kg)
g = 0.0 #acceleration of gravity (m/s^2)
dt = .02 #time step (s)
t_stop = 20.0 #total time (s)
#if the animation is dynamic, set the number of steps based on the stop time and time step
if dyn:
    steps = int(t_stop/dt)

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
d = dfunc.cospath(steps,3,H,.5)
# d1 = dfunc.linepath(steps/3,H*.5)
# d2 = dfunc.pausepath(steps/3,H*.5)
# d3 = dfunc.linepath(steps/3,H,d0=H*.5)
# d=d1+d2+d3



#--------------------check for input errors------------------------------------------------------


#binning related
if bing:
    if i_max >= n_l:
        raise ValueError('The number of bins must be smaller than the number of links')


#damage related
if damage:
    if theta_crit < 0 or theta_crit > 90:
        raise ValueError('The breaking angle must be between 0 and 90 degrees')


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


def Analysis(ea_t, k_t, H, Wm, n_l, d, i_max, theta_crit, m, g, dt, t_stop, dyn, damage, bing, testing):
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
    
    
    #set up stuff for root finding if not using dynamic analysis
    if not dyn:
        #set initial root guesses
        F_old = F_crit #initial guesses
        F_old_avg = F_crit_avg #initial guess for average
        
        #define function to find root of
        def root_g(f,d,L,k,H,ea):
            return L*(1-(4*f)/(k*L)) + (H-L)*math.tanh(f/ea) - d
        
        #define derivative of function
        def root_gprime(f,d,L,k,H,ea):
            return -4/k + ((H-L)/ea)*(1-(math.tanh(f/ea))**2)
    #set up stuff for ODE solving if using dynamic analysis
    else:
        #initialize velocity
        v = [0]*n
        v_avg = 0
        #initialize acceleration
        a = 0
        a_avg = 0
    
    
    #calculate force at each displacement value
    for i in range(len(d)):
        #print i   
        
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
                    
                        #set the guess
                        guess = F_old[j]
                        
                        #set zerotol
                        if L[j] < .3:
                            ztol = .1
                        else:
                            ztol = .01
                        
                        #calculate the force
                        F[j] = MyRootFinding.newton_raphson(root_g,guess,root_gprime,args=(d[i],L[j],k,H,ea),tol=.001,Imax=10*len(d),zerotol=ztol)
                        
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
                    #solve ODE for y with predictor-corrector method
                    
                    #calculate the acceleration
                    #if the link isn't broken
                    if b[j] == 0:
                        #if the argument for atanh steps over the lower boundary, set to lowest possible value to correct
                        if (d[i-1]-y[i-1][j])/(H-L[j]) <= -1:
                            print "invalid atanh argument"
                            a = (k/(4*m))*(L[j]-y[i-1][j]) - (ea/m)*math.atanh(-.9999999999999999) - g
                        #if the argument for atanh steps over the upper boundary, set to highest possible value to correct
                        elif (d[i-1]-y[i-1][j])/(H-L[j]) >= 1:
                            print "invalid atanh argument"
                            a = (k/(4*m))*(L[j]-y[i-1][j]) - (ea/m)*math.atanh(.9999999999999999) - g
                        #otherwise, calculate as normal
                        else:
                            a = (k/(4*m))*(L[j]-y[i-1][j]) - (ea/m)*math.atanh((d[i-1]-y[i-1][j])/(H-L[j])) - g
                    #if the link is broken, the side spring and thus the link provides no force
                    else:
                        #if the argument for atanh steps over the lower boundary, set to lowest possible value to correct
                        if (d[i-1]-y[i-1][j])/(H-L[j]) <= -1:
                            print "invalid atanh argument"
                            a = -(ea/m)*math.atanh(-.9999999999999999) - g
                        #if the argument for atanh steps over the upper boundary, set to highest possible value to correct
                        elif (d[i-1]-y[i-1][j])/(H-L[j]) >= 1:
                            print "invalid atanh argument"
                            a = -(ea/m)*math.atanh(.9999999999999999) - g
                        #otherwise, calculate as normal
                        else:
                            a = -(ea/m)*math.atanh((d[i-1]-y[i-1][j])/(H-L[j])) - g
                     
                    #calculate intermediate v
                    v1 = v[j] + a*dt
                    
                    #calculate intermediate y
                    y1 = y[i-1][j] - v[j]*dt
                    
                    #if y1 has stepped over the boundary for the allowable atanh argument, correct
                    #if y1 is too negative, increase slightly
                    if d[i]-y1 >= H-L[j]:
                        y1 = d[i]-H+L[j]+.001*H
                    #if y1 is too positive, decrease slightly
                    elif d[i]-y1 <= L[j]-H:
                        y1 = d[i]+H-L[j]-.001*H
                    
                    #calculate intermediate acceleration
                    #if the link isn't broken
                    if b[j] == 0:
                        a1 = (k/(4*m))*(L[j]-y1) - (ea/m)*math.atanh((d[i]-y1)/(H-L[j])) - g
                    #if the link is broken (no link force)
                    else:
                        a1 = -(ea/m)*math.atanh((d[i]-y1)/(H-L[j])) - g
                    
                    #apply predictor-corrector formula to find new v
                    v[j] = v[j] + (dt/2)*(a1+a)
                    
                    #apply predictor-corrector formula to find new y
                    y[i][j] = y[i-1][j] - (dt/2)*(v1+v[j])
                    
                    #if y[i][j] has stepped over the boundary for the allowable atanh argument, correct and calculate force
                    #if y[i][j] is too negative, increase slightly
                    if d[i]-y[i][j] >= H-L[j]:
                        y_corr = d[i]-H+L[j]+.001*H
                        #calculate force using the corrected y value
                        F[j] = ea*math.atanh((d[i]-y_corr)/(H-L[j]))
                    #if y[i][j] is too positive, decrease slightly
                    elif d[i]-y[i][j] <= L[j]-H:
                        y_corr = d[i]+H-L[j]-.001*H
                        #calculate force using the corrected y value
                        F[j] = ea*math.atanh((d[i]-y_corr)/(H-L[j]))
                    #if y[i][j] is fine, calculate normally
                    else:
                        #calculate force using the y value
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
                if y[i][j] < 0:
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
                    guess = F_old_avg
                    
                    #calculate the force
                    F_avg = MyRootFinding.newton_raphson(root_g,guess,root_gprime,args=(d[i],L_avg,k,H,ea),tol=.001,Imax=10*len(d),zerotol=.01)
                    
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
                #solve ODE for y with predictor-corrector method
                
                #calculate the acceleration
                #if the link isn't broken
                if b_avg == 0:
                    #if the argument for atanh steps over the lower boundary, set to lowest possible value to correct
                    if (d[i-1]-y_avg[i-1])/(H-L_avg) <= -1:
                        print "invalid atanh argument"
                        a_avg = (k/(4*m))*(L_avg-y_avg[i-1]) - (ea/m)*math.atanh(-.9999999999999999) - g
                    #if the argument for atanh steps over the upper boundary, set to highest possible value to correct
                    elif (d[i-1]-y_avg[i-1])/(H-L_avg) >= 1:
                        print "invalid atanh argument"
                        a_avg = (k/(4*m))*(L_avg-y_avg[i-1]) - (ea/m)*math.atanh(.9999999999999999) - g
                    #otherwise, calculate as normal
                    else:
                        a_avg = (k/(4*m))*(L_avg-y_avg[i-1]) - (ea/m)*math.atanh((d[i-1]-y_avg[i-1])/(H-L_avg)) - g
                #if the link is broken, the side spring and thus the link provides no force
                else:
                    #if the argument for atanh steps over the lower boundary, set to lowest possible value to correct
                    if (d[i-1]-y_avg[i-1])/(H-L_avg) <= -1:
                        print "invalid atanh argument"
                        a_avg = -(ea/m)*math.atanh(-.9999999999999999) - g
                    #if the argument for atanh steps over the upper boundary, set to highest possible value to correct
                    elif (d[i-1]-y_avg[i-1])/(H-L_avg) >= 1:
                        print "invalid atanh argument"
                        a_avg = -(ea/m)*math.atanh(.9999999999999999) - g
                    #otherwise, calculate as normal
                    else:
                        a_avg = -(ea/m)*math.atanh((d[i-1]-y_avg[i-1])/(H-L_avg)) - g
                 
                #calculate intermediate v_avg
                v1_avg = v_avg + a_avg*dt
                
                #calculate intermediate y_avg
                y1_avg = y_avg[i-1] - v_avg*dt
                
                #if y1_avg has stepped over the boundary for the allowable atanh argument, correct
                #if y1_avg is too negative, increase slightly
                if d[i]-y1_avg >= H-L_avg:
                    y1_avg = d[i]-H+L_avg+.001*H
                #if y1_avg is too positive, decrease slightly
                elif d[i]-y1_avg <= L_avg-H:
                    y1_avg = d[i]+H-L_avg-.001*H
                
                #calculate intermediate acceleration
                #if the link isn't broken
                if b_avg == 0:
                    a1_avg = (k/(4*m))*(L_avg-y1_avg) - (ea/m)*math.atanh((d[i]-y1_avg)/(H-L_avg)) - g
                #if the link is broken (no link force)
                else:
                    a1_avg = -(ea/m)*math.atanh((d[i]-y1_avg)/(H-L_avg)) - g
                
                #apply predictor-corrector formula to find new v
                v_avg = v_avg + (dt/2)*(a1_avg+a_avg)
                
                #apply predictor-corrector formula to find new y
                y_avg[i] = y_avg[i-1] - (dt/2)*(v1_avg+v_avg)
                
                #if y_avg[i] has stepped over the boundary for the allowable atanh argument, correct
                #if y_avg[i] is too negative, increase slightly
                if d[i]-y_avg[i] >= H-L_avg:
                    y_avg[i] = d[i]-H+L_avg+.001*H
                #if y_avg[i] is too positive, decrease slightly
                elif d[i]-y_avg[i] <= L_avg-H:
                    y_avg[i] = d[i]+H-L_avg-.001*H
                
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
            if y_avg[i] < 0:
                #set y correctly
                y_avg[i] = 0
                #set indicator to unbuckled
                c_avg = 0
                #recalculate force to reflect unbuckling
                F_avg = ea*math.atanh(d[i]/(H - L_avg))
        
        
        #if the link is completely collapsed
        elif c_avg == 2:#-----------------------------------------------------------
            
            #store y
            y_avg[i] = L_avg
            
            #calculate the force in the link
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
        P[k], P_avg, y, F_crit_avg[k], L = Analysis(ea_t, k_t, H, Wm, n_l, d, i_max, theta_crit, m, g, dt, t_stop, dyn, damage, bing, testing)
#if only running once
else:
    P, P_avg, y, F_crit_avg, L = Analysis(ea_t, k_t, H, Wm, n_l, d, i_max, theta_crit, m, g, dt, t_stop, dyn, damage, bing, testing)
    
    
#plot

#if multirun is on
if multirun:
    for k in range(runs):
        plt.plot(d,P[k],'b',lw=2,alpha=0.2)
        
    plt.axis([0,H+.1*H,0,max(F_crit_avg)*n_l*1.5])
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
        plt.axis([0,H+.1*H,min(P_avg)*1.1,max(P_avg)*1.1])
    plt.xlabel("Displacement (m)")
    plt.ylabel("Force (N)")
    plt.title("Force vs. Displacement")
    plt.legend(["Actual Force Value","Force Value for Average Length"],fontsize = 'medium')



#stop timing code
print "Elapsed time is %f seconds" % (time.clock() - start_time)

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
