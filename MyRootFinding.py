#imports
import math
import matplotlib.pyplot as plt
from numpy import linspace



def newton_raphson(f,x0,fprime,args=(),tol=.001,Imax=1000,zerotol=.01):
#     Uses the Newton Raphson method to find the zero of a function
#     Inputs: function, initial guess, derivative, other arguments (in tuple form), tolerance, and max iterations
#     Outputs: location of zero, if found, or various useful data if zero is not found.

    #initialize
    err = 10*tol #initial error
    I = 0 #iteration counter
    fd0 = fprime(x0,*args)
    adj = 0
    #initialize x
    #if slope at initial guess is very small
    if abs(fd0) < .05:
        #change initial guess
        x = x0 + .5*x0
    else:
        #use given initial guess
        x = x0
            
    while I < Imax and err >= tol:
        #print 'I', I
        #store old value
        x_old = x
        #compute values of function and derivative
        fv = f(x,*args)
        fd = fprime(x,*args)
        #find next guess
        x = x - fv/fd
        
        #calculate error
        err = abs((x - x_old)/x)
        
        
        I += 1 #increment counter
        
        
    
    
    fv = f(x,*args)
    
    if I == Imax:
        if abs(fv)<=zerotol:
            print "root finding: max iterations reached, but abs of function value %f <= zerotol %f" % (fv,zerotol)
            return x
        else:    
            print "root finding: max iterations reached"
            print "root finding: function value of %f at location %f" % (fv,x)
            #plot function and mark initial guess
            y = linspace(0,x0*2,100)
            fpl=[f(z,*args) for z in y]
            print "root finding: max function value: ", max(fpl)
            if abs(max(fpl))<=zerotol:
                print "root finding: function gets closer to zero than +/- %f" %zerotol
                print "root finding: using function value closest to zero"
                fpl_max = max(fpl)
                index = [i for i,j in enumerate(fpl) if j == fpl_max]
                x = y[index]
                print "root finding: function value of %f at %f" % (fpl_max,x)
                return x
            else:
                print "root finding: arguments: ", args
                print "root finding: displaying plot of function"
                plt.plot(y,fpl,[x0,x0],[-.1,.1],'r',[y[0],y[99]],[0,0],'k')
                plt.title("function and initial guess")
                plt.show()
            
    else:
        return x        



def bisection(f,xl,xu,args=(),tol=.001,Imax=1000):
    #this function uses the bisection method to find the root of a function
    #inputs: f, function to find the root of
    #        xl and xu, the lower and upper boundaries of the function
    #        args, optional, the other arguments of f
    #        tol, optional (default .001), the error tolerance
    #        Imax, the maximum number of iterations to run the method
    #outputs: the root of the function
    
    
    #initialize
    #find the value of the function at the lower bound
    fxl = f(xl,*args)
    #set the error
    err = 10*tol
    #initialize counter
    I = 0
    #initialize guess for root
    x = xl
    
    
    #perform bisection method
    while (err > tol and I < Imax):
        #print 'I', I
        #save old guess value
        xold = x
        
        #get new guess value at halfway point
        x = (xl + xu)/2.
        #print 'x', x
        #evaluate function at root
        fx = f(x,*args)
        
        #check to see whether the guess has the same sign as the lower bound
        fsign = fxl*fx
        
        #if both have the same sign
        if fsign > 0:
            #move lower bound to guess value
            xl = x
            #save new lower bound function value
            fxl = fx
        #if both have different signs
        elif fsign < 0:
            #move upper bound to guess value
            xu = x
        #if the guess is actually the root
        else:
            #stop the loop
            break
        
        #increment counter
        I += 1
        
        #calculate error
        err = abs((x-xold)/x)
        
        
    if I == Imax:
        print 'root not found'
    else:
        return x
            
    
    
    
    
    
    
    
    
    