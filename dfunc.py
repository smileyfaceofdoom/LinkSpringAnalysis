'''Katharin Jensen
These functions create a few displacement paths for the link-spring system'''

import numpy as np
import matplotlib.pyplot as plt

#basic straight line
def linepath(s,H,d0=0):
    #This function creates a steady straight line displacement path
    #input: The number of values (s), the total displacement (H), starting displacement (d0) (default 0)
    #output: list of displacement values (d)
    
    d = np.linspace(d0,H,s)
    d=list(d)
    
    return d


#v shaped paths
def vpath(s,vs,H,per,d0=0):
    #This function creates v-shaped displacement paths, like a sawtooth.
    #inputs: Number of values (s), number of v's (vs)n (can be whole numbers e.g. 1 or half numbers e.g. 2.5), 
    #        maximum possible displacement (H), fraction of H to displace to (per), starting displacement (d0) (default 0)
    #outputs: list of displacements (d)
    
    #create list of step values
    s_vals = range(s)
    #find what fraction of s_vals goes into each line segment
    s_frac = int(np.floor(s/(2*vs)))
    #initialize d
    d = [0]*(s-s%s_frac)
    #find slope
    m = (H*per-d0)/((s_vals[len(s_vals)-1]-s%s_frac)/(2*vs))
    #find equation for each line segment, add to d
    for k in range(int(2*vs)):
        #find x- and y-intercepts
        #for positive slope segments
        if k%2 == 0:
            xint = s_frac*k
            yint = -m*xint
            d_seg = [m*x + yint + d0 for x in s_vals[k*s_frac:(k+1)*s_frac]]
        #for negative slope segments    
        else:
            xint = s_frac*(k+1)
            yint = m*xint
            d_seg = [-m*x + yint + d0 for x in s_vals[k*s_frac:(k+1)*s_frac]]
        #add d_seg to d
        d[k*s_frac:(k+1)*s_frac] = d_seg
    d=list(d)
        
    return d


#sin/cos shaped paths
def cospath(s,w,H,per,d0=0):
    #This function makes a cosine wave shaped path
    #inputs: number of values (s), number of waves (w), 
    #        maximum possible displacement (H), fraction of H to displace t0 (per), starting displacement (d0) (default 0)
    #outputs: list of displacement values (d)
    
    #make sure w is a float
    w = float(w)
    #make s list
    s_vals = range(s)
    #find displacement values
    d = [-((H*per-d0)/2)*(np.cos((float(sv)/s)*w*2*np.pi) - 1) + d0 for sv in s_vals]
    d=list(d)
    
    return d


#put a pause in the path
def pausepath(s,dval):
    #this function can be used in combination with other functions to put a pause in the displacement
    #inputs: number of values (s), displacement value to pause at (dval)
    #outputs: list of displacement values (d)
    
    #make d list
    d = [dval for sval in range(s)]
    
    return d
    

#tests

#s = range(999/3)
#d = linepath(int(20/(3*.02)),1.0,d0=1.0/2)

# d = vpath(1000,1.5,1.0,.5,d0=.25)
# s = range(len(d))

# s = range(1000)
# d = cospath(1000,2.3,1.0,.5,.25)

# s = range(1000)
# d = pausepath(1000,.75)
  
# plt.plot(s,d)
# plt.show()



