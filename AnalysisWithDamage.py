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
ea_t = 10.0 #top spring constant (N)
k_t = 100.0 #side spring constant (N/m)
H = 1.0 #Total height (m)
m = 4.0 #weibull modulus
n_l = 3 #number of links
dmax = .5 #displacement before unloading (Set to H if unloading isn't desired)
theta_crit = 50 #angle from vertical at which the links "break", degrees (set to 180 if breaking isn't desired)
i_max = 21 #number of nodes for binning
bing = False #apply binning? True or False

#makes sure the number of links to break is not larger than the total number of links
if dmax > H:
    raise ValueError('The displacement before unloading can\'t be larger than the maximum displacement')
#make sure theta_crit isn't negative or a value larger than the largest possible rotation angle (except the value for not breaking)
if (theta_crit > 90 and theta_crit != 180) or theta_crit < 0:
    raise ValueError('Breaking angle should be either positive and <= 90 degrees or 180 degrees for no breaking.')

#normalize
ea = ea_t/n_l
k = k_t/n_l

#convert theta to radians
theta_crit = theta_crit*math.pi/180

#generate list of random numbers
R = [random() for x in range(n_l)]
#for testing purposes
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

#initialize displacement list
d = [0]
ldng = 1 #set whether loading (1) or unloading (0)
#initialize y-values2
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
    
    
    #for the average length
    #if the links haven't broken
    if b_avg == 0:
        if c_avg == 0 and d[i] >= (H-L_avg):
                #if so, set indicator to buckling
                c_avg = 1
            
            #if it hasn't buckled
        if c_avg == 0:
            #calculate force in link
            F_avg = ea*math.atanh(d[i]/(H-L_avg))
            #store y value
            y_avg.append(0)
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
            y_avg.append(L_avg - 4*F_avg/k)
            if y_avg[i] > (L_avg - 4*F_avg/k):
                y_avg[i] = L_avg
            #check if link has broken
            if (L_avg-y_avg[i])/L_avg <=1 and (L_avg-y_avg[i])/L_avg >= -1:
                theta = math.acos((L_avg-y_avg[i])/L_avg)
            else:
                theta = math.pi/2    
            #check if loading
            if ldng == 1:
                if theta > theta_crit:
                    b_avg = 1
                    #check if top spring can extend all the way
                    if d[i]>L_avg:
                        c_avg = 2
            else:#if unloading
                #check if the links have unbent again
                    if y_avg[i] < 0:
                        y_avg[i] = 0
                        #set force correctly
                        F[j] = ea*math.atanh(d[i]/(H-L_avg))
                        c_avg = 0
                            
                    
            #check if link has completely collapsed
            if  y_avg[i] >= L_avg:
                c_avg = 2 #set indicator
                F_avg = ea*math.atanh((d[i] - L_avg)/(H - L_avg))
        #if it's completely collapsed
        else:
            #calculate link force
            F_avg = ea*math.atanh((d[i] - L_avg)/(H - L_avg))
            #store y
            y_avg.append(L_avg)
    else: #if links have broken
        #if top spring is not fully extended
            if c[j] == 2:
                #calculate link force
                F_avg = ea*math.atanh((d[i] - L_avg)/(H - L_avg))
                #store y
                y_avg.append(L_avg)
                
                #see if spring is fully extended
                if d[i] <= L_avg:
                    #set indicator
                    c_avg = 1
            elif c_avg == 1: #if top spring is fully extended
                #force is zero
                F_avg = 0
                #store y value
                y_avg.append(d[i])
            else:
                print "springs are extending"        
    
    #find P_avg
    P_avg.append(F_avg*n_l)
    
    
    
    #increment counter
    i += 1
    

#plot
plt.plot(d,P,d,P_avg,'r--')
plt.axis([0,H+.1*H,0,F_crit_avg*n_l*1.1])
plt.xlabel("Displacement (m)")
plt.ylabel("Force (N)")
plt.title("Force vs. Displacement")
plt.legend(["Actual Force Value","Force Value for Average Length"],fontsize = 'medium')

print "Elapsed time is %f seconds" % (time.clock() - start_time)
plt.show()

#animate: note, binning must be turned off for animations.
#Would not recommend using for more than 5 links, as it gets really small after that.

#initialization
#initialize figure
fig = plt.figure()
#set axes
ax1 = fig.add_subplot(211)
ax1.set_xlim(-H*.6,H*n+.6*H)
ax1.set_ylim(-.1*H,H*(n)/2+.1*H)
ax2 = fig.add_subplot(212)
ax2.set_xlim(0,H+.1*H)
ax2.set_ylim(-.1,F_crit_avg*n_l*1.1)
#initialize patch vertices and move codes
verts = [[-10.,-10.],[-10.,-9.],[-9.,-9.],[-9.,-10.],[-10.,-10.]]
codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]
#set path
path = Path(verts,codes)
#create patch
patch = patches.PathPatch(path,facecolor='white',edgecolor='white')
#create screen clearer
clr = patches.Rectangle((-H*.6,-.1*H),H*n+.6*H+H*.6,H*(n)/2+.1*H+.1*H,facecolor='white',edgecolor='white')
#draw initial (blank) lines and patches
lineb, = ax1.plot([],[],'k',lw=7)
linesp, = ax1.plot([],[],'k',lw=2)
dots, = ax1.plot([],[],'k.')
ax1.add_patch(patch)
linePd, = ax2.plot([],[])
linePd_avg, = ax2.plot([],[],'r--')
ax2.legend(["Actual Force Value","Force Value for Average Length"],fontsize='small')
ax1.set_title('Link Spring System')
ax2.set_title('Force vs. Displacement')
ax2.set_xlabel('Displacement (m)')
ax1.set_ylabel('m')
ax2.set_ylabel('Force (N)')

#set initialization function for animation
def init():
    #clear patch
    ax1.add_patch(clr)
    #empty lines
    linesp.set_data([],[])
    lineb.set_data([],[])
    dots.set_data([],[])
    linePd.set_data([],[])
    linePd_avg.set_data([],[])
    #blank patch
    ax1.add_patch(patch)
    #return data
    return clr,lineb,patch,linesp,dots,linePd_avg,linePd

def animate(i):
    #clear frame
    clr = patches.Rectangle((-H*.6,-.1*H),H*n+.6*H+H*.6,H*(n)/2+.1*H+.1*H,facecolor='white',edgecolor='white')
    ax1.add_patch(clr)
    #get frame data from ani function
    path,xybar,spr1x,spr1y,dot = ani.ani(L,y[i],H,n,d[i])
    #set top bar data
    lineb.set_data(xybar[0],xybar[1])
    #set spring data
    linesp.set_data(spr1x,spr1y)
    #set patches for links
    patch = patches.PathPatch(path,facecolor='gray',lw=2)
    ax1.add_patch(patch)
    #add dots
    dots.set_data(dot[0],dot[1])
    #plot force vs. displacement
    linePd.set_data(d[:i],P[:i])
    #plot average force vs. displacement
    linePd_avg.set_data(d[:i],P_avg[:i])
    return clr,lineb,patch,linesp,dots,linePd_avg,linePd

anim = animation.FuncAnimation(fig,animate,init_func=init,frames=(len(d)-1),interval=20,blit=True,repeat=False)
plt.show()

start_ani = time.clock()
print "begin saving animation"
anim.save('3linkwithbreaking.gif', writer='imagemagick')
print "done saving animation"
print "animation save time is %f seconds" % (time.clock() - start_ani)
