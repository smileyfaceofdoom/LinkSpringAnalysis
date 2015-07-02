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
from RK4Dyn import RK4
import dfunc

#set parameters
ea_t = 10.0 #top spring constant thing
k_t = 100.0 #side spring constant
H = 1.0 #max height
Wm = 4.0 #Weibull modulus
m = 1.0 #mass (kg)
g = 0 #acceleration of gravity
n_l = 5 #number of links
dt = .02 #time step (s)
t_stop = 20.0 #time to stop (s)
xspacing = .5*H #spacing for animation
testing = True
auto_u_fix = True #if the displacement goes higher than the allowed, set to max allowed value? Yes if True. If False, will break out of loop.
save_ani = False


#create list of displacements
u = dfunc.cospath(int(t_stop/dt),3,H,.5)


#normalize
ea = ea_t/n_l
k = k_t/n_l


#generate random numbers
if testing:
    #for testing purposes
    R = [0.1406236293192208, 0.557452455836349, 0.4018884612118805, 0.8610494090625574, 0.005928894753714831]
    #R = [.5]
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
print "L", L
#L = [.7]
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
print "Fcrit", F_crit

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
            #print F[j]
            #check to see if link has buckled
            if (F[j] + g*m) >= F_crit[j]:
                c[j] = 1 #set indicator
        #if the link is currently buckling:
        elif c[j] == 1:
            #solve ODE for y with predictor-corrector method
            if (u[i-1]-y[i-1][j])/(H-L[j]) <= -1:
                print "invalid atanh argument"
                a = (k/(4*m))*(L[j]-y[i-1][j]) - (ea/m)*math.atanh(-.9999999999999999) - g
                print "a", a
            elif (u[i-1]-y[i-1][j])/(H-L[j]) >= 1:
                print "invalid atanh argument"
                a = (k/(4*m))*(L[j]-y[i-1][j]) - (ea/m)*math.atanh(.9999999999999999) - g
                print "a", a
            else:
                a = (k/(4*m))*(L[j]-y[i-1][j]) - (ea/m)*math.atanh((u[i-1]-y[i-1][j])/(H-L[j])) - g
                print "a", a
            v1 = v + a*dt
            y1 = y[i-1][j] - v*dt
            #print "y1", y1
            
            if u[i]-y1 >= H-L[j]:
                y1 = u[i]-H+L[j]+.001
                #print "adj y1", y1
            elif u[i]-y1 <= L[j]-H:
                y1 = u[i]+H-L[j]-.001
                #print "adj y1", y1
            
            a1 = (k/(4*m))*(L[j]-y1) - (ea/m)*math.atanh((u[i]-y1)/(H-L[j])) - g
            v = v + (dt/2)*(a1+a)
            y[i][j] = y[i-1][j] - (dt/2)*(v1+v)
            
            #y[i][j],v1 = RK4(u[i-1],u[i],y[i-1][j],v,H,L[j],k,ea,m)
            
            #set new v
            #v = v1
            
            if u[i]-y[i][j] >= H-L[j]:
                y[i][j] = u[i]-H+L[j]+.001
                #print "adj y[i][j]", y[i][j]
            elif u[i]-y[i][j] <= L[j]-H:
                y[i][j] = u[i]+H-L[j]-.001
                #print "adj y[i][j]", y[i][j]
            
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
#             print "c2 u", u[i]
#             print "c2 L", L[j]
#             print "atanh", (u[i] - L[j])/(H - L[j])
            
            #store y
            y[i][j] = L[j]
            #set velocity back to zero
            v = 0
            
            if (u[i] - L[j])/(H - L[j]) <= -1:
                F[j] = ea*math.atanh(-.9999999999999999)
            else:
                F[j] = ea*math.atanh((u[i] - L[j])/(H - L[j]))
            #check to see if link is uncollapsed
            if F[j] <= -g*m:
                c[j]=1
                #F[j] = F[j-1]
            #else:
                
    #calculate total force
    P[i] = sum(F)

plt.plot(u,P)
plt.axis([0,H+.1*H,0,16])
plt.xlabel("Displacement (m)")
plt.ylabel("Force (N)")
plt.title("Force vs. Displacement")
plt.show()


P_avg = [0]*len(u)
anim = ani.ani(L,y,H,n,u,P,P_avg,xspacing,F_crit_avg,save_ani,frameskip=1,dyn=True,t_step = dt)

if save_ani:
    #save animation as GIF
    #IMPORTANT NOTE: This will probably take 1-2 hours. Plan accordingly.
    start_ani = time.clock()
    print "begin saving animation"
    anim.save('5linkdynex3.gif', writer='imagemagick')
    print "done saving animation"
    print "animation save time is %f seconds" % (time.clock() - start_ani)



