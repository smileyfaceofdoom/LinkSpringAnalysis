"""Katharin Jensen
This code will analyze the motion of a dynamic link spring system under a load"""


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
import ani

#set parameters
ea_t = 10.0 #top spring constant thing
k_t = 100.0 #side spring constant
H = 1.0 #max height
Wm = 4.0 #Weibull modulus
m = 1 #mass (kg)
n_l = 2 #number of links
dt = .02 #time step (s)
t_stop = 20.0 #time to stop (s)
xspacing = .5*H #spacing for animation
testing = False
auto_u_fix = True #if the displacement goes higher than the allowed, set to max allowed value? Yes if True. If False, will break out of loop.

#set displacement function u(t) (either a function with input t or a list of displacements at even time intervals)
def ufunc(t):
    
    u = math.sin(t*6*math.pi/20.0 - math.pi/2)/4 + .25
    
    return u

#create list of times
t_list = linspace(0,t_stop,(t_stop/dt))


#create list of displacements
u = [ufunc(t) for t in t_list]
print u
#normalize
ea = ea_t/n_l
k = k_t/n_l


#generate random numbers
if testing:
    #for testing purposes
    R = [.4,.3]
else:
    #generate list of random numbers
    R = [random() for x in range(n_l)]
    
    
#create length distribution
def length_distribution(r):
    l = [H*((math.log(1/num)/math.log(2))**(1.0/Wm)) for num in r]
    L = [H*(1-math.exp(-num)) for num in l]
    return L

   
#find length distribution
L = length_distribution(R)
#find average length
L_avg = mean(L)
n = n_l    
    
    
#initialize y list
y = [[0]*n for x in range(len(u))]
y_avg = [0]*(len(u))
#initialize mass velocity
v = 0
#initialize mass acceleration
a = 0
#initialize force list
P = [0]*(len(u))
P_avg = [0]*(len(u))
#set indicator for whether links have buckled: 0 = no, 1 = in progress, 2 = yes
c = [0]*n
c_avg = 0  

#find critical buckling values
F_crit = [k*Li/4 for Li in L]
F_crit_avg = k*L_avg/4


#calculate the force at each displacement value
for i in range(len(u)):
    #make sure the displacement is not larger than it's allowed to be
    if u[i] >= H:
        if auto_u_fix:
            print "Error: displacement value larger than max allowed."
            print "Setting displacement to maximum allowed value"
            u[i] = .99*H
        else:
            print "Error: displacement value larger than max allowed."
            break
    #initialize force list
    F = [0]*n
    
    #for each link
    for j in range(n):
        #if the link hasn't buckled
        #check to see if displacement is past top of link
        if c[j] == 0 and u[i] >= (H-L[j]):
            #if so, set indicator to buckling
            c[j] = 1
            
        #if link still hasn't buckled
        if c[j] == 0:
            #calculate force in link
            F[j] = ea*math.atanh(u[i]/(H-L[j]))
            y[i][j] = 0
            #check to see if link has buckled
            if (F[j] + 9.81*m) >= F_crit[j]:
                c[j] = 1 #set indicator
        #if the link is currently buckling:
        elif c[j] == 1:
            #solve ODE for y with predictor-corrector method
            if (u[i-1]-y[i-1][j])/(H-L[j]) <= -1:
                a = (k/(4*m))*(L[j]-y[i-1][j]) - (ea/m)*math.atanh(-.9999999999999999) - 9.81
            else:
                a = (k/(4*m))*(L[j]-y[i-1][j]) - (ea/m)*math.atanh((u[i-1]-y[i-1][j])/(H-L[j])) - 9.81
            v1 = v + a*dt
            y1 = y[i-1][j] - v1*dt
            print "y1", y1
            
            if u[i]-y1 >= H-L[j]:
                y1 = u[i]-H+L[j]+.001
                print "adj y1", y1
            elif u[i]-y1 <= L[j]-H:
                y1 = u[i]+H-L[j]-.001
                print "adj y1", y1
            print "u", u[i]
            print "L", L[j]
            print "atanh", (u[i]-y1)/(H-L[j])
            
            a1 = (k/(4*m))*(L[j]-y1) - (ea/m)*math.atanh((u[i]-y1)/(H-L[j])) - 9.81
            v = v + (dt/2)*(a1+a)
            y[i][j] = y[i-1][j] - (dt/2)*(v1+v)
            
            if u[i]-y[i][j] >= H-L[j]:
                y[i][j] = u[i]-H+L[j]+.001
                print "adj y[i][j]", y[i][j]
            elif u[i]-y[i][j] <= L[j]-H:
                y[i][j] = u[i]+H-L[j]-.001
                print "adj y[i][j]", y[i][j]
            
            #calculate force
            F[j] = ea*math.atanh((u[i]-y[i][j])/(H-L[j]))
            
            #check if link is still buckling
            #if it's collapsed
            if y[i][j] >= L[j]:
                #calculate force
                F[j] = ea*math.atanh((u[i] - L[j])/(H - L[j]))
                #set y value
                y[i][j] = L[j]
                #set indicator
                c[j] = 2
            #if it's gone back to unbuckled
            elif y[i][j] <= 0:
                #calculate force
                F[j] = ea*math.atanh(u[i]/(H-L[j]))
                #set y value
                y[i][j] = 0
                #set indicator
                c[j] = 0
        #if link is completely collapsed
        else:
            #calculate link force
            print "c2 u", u[i]
            print "c2 L", L[j]
            print "atanh", (u[i] - L[j])/(H - L[j])
            
            #store y
            y[i][j] = L[j]
            #set velocity back to zero
            v = 0
            #check to see if link is uncollapsed
            if u[i]-L[j] <= -(H-L[j]):
                c[j]=1
                F[j] = F[j-1]
            else:
                F[j] = ea*math.atanh((u[i] - L[j])/(H - L[j]))
    #calculate total force
    P[i] = sum(F)

plt.plot(u,P)
plt.axis([0,H+.1*H,-30,30])
plt.xlabel("Displacement (m)")
plt.ylabel("Force (N)")
plt.title("Force vs. Displacement")
plt.show()

save_ani = False
P_avg = [0]*len(u)
anim = ani.ani(L,y,H,n,u,P,P_avg,xspacing,F_crit_avg,save_ani,dyn=True)





