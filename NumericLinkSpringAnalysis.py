""" Katharin Jensen
This code solves for the force-displacement curve for a system of links and springs (described elsewhere).
It then creates a force vs. displacement plot."""

#imports
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

#begin timing code
start_time = time.clock()


#initialize
ea_t = 10.0 #top spring constant
k_t = 100.0 #side spring constant
H = 1.0 #Total height
m = 4.0 #weibull modulus
n_l = 5 #number of links
i_max = 6 #number of nodes for binning
bing = False #apply binning? True or False
make_ani = True #create animation? True or False
save_ani = False #save animation as .gif? True or False
testing = True #Use a non-random list of R values for testing purposes
#set link spacing for animation
xspacing = .5*H
    
#normalize
ea = ea_t/n_l
k = k_t/n_l

#generate list of random numbers
if testing:
    #for testing purposes
    #R = linspace(.3,.8,n_l)
    #pseudo-random 5-link
    R = [0.1406236293192208, 0.557452455836349, 0.4018884612118805, 0.8610494090625574, 0.005928894753714831]
else:
    #generate list of random numbers
    R = [random() for x in range(n_l)]


#create length distribution
def length_distribution(r):
    l = [H*((math.log(1/num)/math.log(2))**(1.0/m)) for num in r]
    L = [H*(1-math.exp(-num)) for num in l]
    return L

#if binning is applied
if bing:
    #apply binning
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
    #L=[.5]
    #find average length
    L_avg = mean(L)
    n = n_l


#setup for force calculations

#make displacement distribution
if save_ani:
    d = linspace(0,H,200)
else:
    d = linspace(0,H,1000)
#initialize y-values
y = [[0]*n for x in range(len(d)-1)]
y_avg = [0]*(len(d)-1)
#initialize total force list
P = [0]*(len(d)-1)
P_avg = [0]*(len(d)-1)
#set indicator for whether links have buckled: 0 = no, 1 = in progress, 2 = yes
c = [0]*n
c_avg = 0

#find critical buckling values
F_crit = [k*Li/4 for Li in L]
F_crit_avg = k*L_avg/4

#define stuff for the newton function
F_old = F_crit #initial guesses
F_old_avg = F_crit_avg #initial guess for average
#define function
def g(f,d,L,k,H,ea):
    return L*(1-(4*f)/(k*L)) + (H-L)*math.tanh(f/ea) - d
#define derivative of function
def gprime(f,d,L,k,H,ea):
    return -4/k + ((H-L)/ea)*(1-(math.tanh(f/ea))**2)


#calculate force at each displacement value (except the last, because that one's infinity)
for i in range(len(d)-1):
    #initialize force list
    F = [0]*n
    #for each link
    for j in range(n):
        #if the link hasn't buckled:
        #check to see if displacement is past top of link
        if c[j] == 0 and d[i] >= (H-L[j]):
            #if so, set indicator to buckling
            c[j] = 1
        
        #if it hasn't buckled
        if c[j] == 0:
            #calculate force in link
            F[j] = ea*math.atanh(d[i]/(H-L[j]))
            y[i][j] = 0
            #check to see if link has buckled
            if F[j] >= F_crit[j]:
                c[j] = 1 #set indicator
        #if it's currently buckling
        elif c[j] == 1: 
            #numerically calculate link force
            guess = F_old[j]
                      
            F[j] = MyRootFinding.newton_raphson(g,guess,gprime,args=(d[i],L[j],k,H,ea),tol=.001,Imax=10*len(d),zerotol=.01)
            #assign current value as next guess value
            F_old[j] = F[j]         
            #find and store y value
            y[i][j] = L[j] - 4*F[j]/k
            #check if link has completely collapsed
            if y[i][j] >= L[j]:
                c[j] = 2 #set indicator
                F[j] = ea*math.atanh((d[i] - L[j])/(H - L[j]))
        #if it's completely collapsed
        else:
            #calculate link force
            F[j] = ea*math.atanh((d[i] - L[j])/(H - L[j]))
            #store y
            y[i][j] = L[j]
    
    #sum forces to find total force
    #if binning is applied
    if bing:
        P[i] = 0 #initialize force sum
        for j in range(n):
            P[i] += F[j]*W_b[j]*n_l #sum weighted forces
    #if binning is not applied        
    else:    
        P[i] = sum(F)    
    
    
    #for the average length
    if c_avg == 0 and d[i] >= (H-L_avg):
            #if so, set indicator to buckling
            c_avg = 1
        
        #if it hasn't buckled
    if c_avg == 0:
        #calculate force in link
        F_avg = ea*math.atanh(d[i]/(H-L_avg))
        #check to see if link has buckled
        if F_avg >= F_crit_avg:
            c_avg = 1 #set indicator
    #if it's currently buckling
    elif c_avg == 1:
        #numerically calculate link force
        F_avg = MyRootFinding.newton_raphson(g,F_old_avg,gprime,args=(d[i],L_avg,k,H,ea),tol=1.0e-4,Imax=5*len(d),zerotol=.01)
        #assign current value as next guess value
        F_old_avg = F_avg         
        #find and store y value
        y_avg[i] = L_avg - 4*F_avg/k
        #check if link has completely collapsed
        if  y_avg[i] >= L_avg:
            c_avg = 2 #set indicator
            F_avg = ea*math.atanh((d[i] - L_avg)/(H - L_avg))
    #if it's completely collapsed
    else:
        #calculate link force
        F_avg = ea*math.atanh((d[i] - L_avg)/(H - L_avg))
        #store y
        y_avg[i] = L_avg
    
    #find P_avg
    P_avg[i] = F_avg*n_l
    

#plot
plt.plot(d[:(len(d)-1)],P,d[:(len(d)-1)],P_avg,'r--')
plt.axis([0,H+.1*H,0,F_crit_avg*n_l*1.1])
plt.xlabel("Displacement (m)")
plt.ylabel("Force (N)")
plt.title("Force vs. Displacement")
plt.legend(["Actual Force Value","Force Value for Average Length"])

print "Elapsed time is %f seconds" % (time.clock() - start_time)
plt.show()

if make_ani:
    #animate: note, binning must be turned off for animations.
    #Would not recommend using for more than 5 links, as it gets really small after that.
    #perform figure setup
    anim = ani.ani(L,y,H,n,d,P,P_avg,xspacing,F_crit_avg,save_ani,fps=50,frameskip=1,dyn=False,t_step = .1)

if save_ani:
    #save animation as GIF
    #IMPORTANT NOTE: This will probably take 1-2 hours. Plan accordingly.
    start_ani = time.clock()
    print "begin saving animation"
    anim.save('5linkex2.gif', writer='imagemagick')
    print "done saving animation"
    print "animation save time is %f seconds" % (time.clock() - start_ani)
