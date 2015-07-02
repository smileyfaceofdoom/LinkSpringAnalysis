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
n_l = 10 #number of links
runs = 100 #number of times to run
dmax = 1
theta_crit = 50
i_max = 21 #number of nodes for binning
bing = False #apply binning? True or False
make_ani = False #create animation? True or False
save_ani = False #save animation as .gif? True or False
testing = False #Use a non-random list of R values for testing purposes
#set link spacing for animation
xspacing = .5*H
 
def NLSA(ea_t,k_t,H,m,n_l):   
    #normalize
    ea = ea_t/n_l
    k = k_t/n_l
    
    #generate list of random numbers
    if testing:
        #for testing purposes
        R = linspace(.3,.8,n_l)
        #pseudo-random 5-link
        #R = [0.1406236293192208, 0.557452455836349, 0.4018884612118805, 0.8610494090625574, 0.005928894753714831]
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
        
    return d,P,F_crit_avg


def NLSADmg(ea_t,k_t,H,m,n_l,dmax,theta_crit):
    #normalize
    ea = ea_t/n_l
    k = k_t/n_l
    
    #convert theta to radians
    theta_crit = theta_crit*math.pi/180
    
    #generate list of random numbers
    R = [random() for x in range(n_l)]
    #for testing purposes
    #R = linspace(.3,.9,n_l)
    #pseudo-random 5-link
    #R = [0.1406236293192208, 0.557452455836349, 0.4018884612118805, 0.8610494090625574, 0.005928894753714831]
    
    
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
    
    #initialize displacement list
    d = [0]
    ldng = 1 #set whether loading (1) or unloading (0)
    #initialize y-values
    y = [[0]*n]
    y_avg = [0]
    #initialize total force list
    P = [0]
    P_avg = [0]
    
    #set indicator for whether links have buckled broken: 0 = no, 1 = buckling, 2 = completely collapsed
    c = [0]*n
    c_avg = 0
    #set indicator for whether links have broken: 0 = no, 1 = yes
    b = [0]*n
    b_avg = 0
    
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
    i = 1 #set counter
    #loop will continue until d comes back up to its starting position, or until it goes through 100000 iterations
    #(to make sure it'll stop even if it screws up)
    while i <= 100000:
        #set new d value (incrementing by thousandths of H)
        #if not enough links have broken yet:
        if d[i-1] >= dmax:
            ldng = 0
        if save_ani:
            if ldng == 1: #if loading    
                #increment d forward
                d.append(i*.005*H)
            else: #if unloading
                #increment d backward
                d.append(d[i-1]-.005*H)
        else:
            if ldng == 1: #if loading    
                #increment d forward
                d.append(i*.001*H)
            else: #if unloading
                #increment d backward
                d.append(d[i-1]-.001*H)
        #stop if d has displaced further than H
        if d[i] >= H:
            #make d the same length as the force lists
            d.pop()
            break
        #stop if d becomes negative
        if d[i] < 0:
            #make d the same length as the force lists
            d.pop()
            break
        #add next list of y-values
        y.append([])        
        #initialize force list
        F = [0]*n
        #for each link
        for j in range(n):
            #if link hasn't broken:
            if b[j] == 0:
                #if the link hasn't buckled:
                #check to see if displacement is past top of link
                if c[j] == 0 and d[i] >= (H-L[j]):
                    #if so, set indicator to buckling
                    c[j] = 1
                
                #if it hasn't buckled
                if c[j] == 0:
                    #calculate force in link
                    F[j] = ea*math.atanh(d[i]/(H-L[j]))
                    y[i].append(0)
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
                    y[i].append(L[j] - 4*F[j]/k)
                    if y[i][j] > L[j]:
                        y[i][j] = L[j]
                    #if still loading
                    #(if not, d is decreasing, and the remaining links aren't going to break and don't need to be checked)
                    if ldng == 1:
                        #check if link has broken
                        #small links will have very large jumps in y, so only calculate acos if theta < 90 degrees
                        if (L[j]-y[i][j])/L[j] <= 1 and (L[j]-y[i][j])/L[j] >= -1:
                            theta = math.acos((L[j]-y[i][j])/L[j])
                        else: #otherwise theta should be 90 degrees
                            theta = math.pi/2 
                               
                        if theta >= theta_crit:
                            b[j] = 1
                            #check whether top spring has enough room to extend fully
                            if d[i]>L[j]:
                                c[j] = 2
                        #check if link has completely collapsed
                        if y[i][j] >= L[j]:
                            c[j] = 2 #set indicator
                            F[j] = ea*math.atanh((d[i] - L[j])/(H - L[j]))       
                    else: #if unloading
                        #check if the links have unbent again
                        if y[i][j] < 0:
                            y[i][j] = 0
                            #set force correctly
                            F[j] = ea*math.atanh(d[i]/(H-L[j]))
                            c[j] = 0        
                #if it's completely collapsed
                else:
                    #calculate link force
                    F[j] = ea*math.atanh((d[i] - L[j])/(H - L[j]))
                    #store y
                    y[i].append(L[j])
            else: #if link has broken
                #if top spring is not fully extended
                if c[j] == 2:
                    #calculate link force
                    F[j] = ea*math.atanh((d[i] - L[j])/(H - L[j]))
                    #store y
                    y[i].append(L[j])
                    #check if spring has fully extended
                    if d[i]<=L[j]:
                        #set indicator
                        c[j] = 1
                elif c[j] == 1: #if top spring is fully extended
                    #force is zero
                    F[j] = 0
                    #store y value
                    y[i].append(d[i])
                    #see if spring is being compressed
                    if d[i] > L[j]:
                        #set indicator
                        c[j] = 2
                else:
                    print "spring %d is extending" % j         
        
    
        #sum forces to find total force
        #if binning is applied
        if bing:
            Psum = 0 #initialize force sum
            for j in range(n):
                Psum += F[j]*W_b[j]*n_l #sum weighted forces
            P.append(Psum)    
        #if binning is not applied        
        else:    
            P.append(sum(F))
            
        #increment counter
        i += 1    
        
    return d,P,F_crit_avg    
    
    
    
d = [[] for i in range(runs)]   
P = [[] for i in range(runs)]
F_crit_avg = [0]*runs
for k in range(runs):
    print k
    d[k],P[k],F_crit_avg[k] = NLSA(ea_t,k_t,H,m,n_l)#,dmax,theta_crit)
    d[k] = list(d[k])
    d[k].pop()
    

#plot

for k in range(runs):
    plt.plot(d[k],P[k],'b',lw=2,alpha=0.2)
    
plt.axis([0,H+.1*H,0,max(F_crit_avg)*n_l*1.5])
plt.xlabel("Displacement (m)")
plt.ylabel("Force (N)")
plt.title("%d %d-link Force vs. Displacement Curves" % (runs, n_l))
#plt.legend(["%d links" % n_l[0],"%d links" % n_l[1],"%d links" % n_l[2],"%d links" % n_l[3],"%d links" % n_l[4],"%d links" % n_l[5]])

print "Elapsed time is %f seconds" % (time.clock() - start_time)
plt.show()

if make_ani:
    #animate: note, binning must be turned off for animations.
    #Would not recommend using for more than 5 links, as it gets really small after that.
    #perform figure setup
    anim = ani.ani(L,y,H,n,d,P,P_avg,xspacing,F_crit_avg,save_ani)

if save_ani:
    #save animation as GIF
    #IMPORTANT NOTE: This will probably take 1-2 hours. Plan accordingly.
    start_ani = time.clock()
    print "begin saving animation"
    anim.save('5linkex2.gif', writer='imagemagick')
    print "done saving animation"
    print "animation save time is %f seconds" % (time.clock() - start_ani)
