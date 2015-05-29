
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
import MyRootFinding
import binning

#begin timing code
start_time = time.clock()


#initialize
ea_t = 10.0 #top spring constant
k_t = 100.0 #side spring constant
H = 1.0 #Total height
m = 4.0 #weibull modulus
n_l = 300 #number of links
i_max = 51 #number of nodes for binning
bing = False #apply binning? True or False

#normalize
ea = ea_t/n_l
k = k_t/n_l


#generate list of random numbers
R = [random() for x in range(n_l)]
#for testing purposes, to compare to MATLAB code
#R = linspace(.3,.9,n_l)


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
    #find average length
    L_avg = mean(L)
    n = n_l


#setup for force calculations

#make displacement distribution
d = linspace(0,H,1000)
#initialize y-values
y = [[0]*n]*(len(d)-1)
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


