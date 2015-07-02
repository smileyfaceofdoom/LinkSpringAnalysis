"""This will animate the link spring system, in conjunction with matplotlib.animation"""
#imports
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.animation as animation
import numpy as np
    
def ani(L,y,H,n,d,P,P_avg,xspacing,F_crit_avg,save_ani,fps=50,frameskip=1,dyn=False,t_step = .1):
    #this function animates the link-spring system.
    #inputs: list of lengths, ALL y values, H, number of links, ALL d values, ALL force values, ALL average force values, spacing between links, F_crit_avg
        
    def draw(L,y,H,n,d,xspacing):
        #This function calculates all the numbers necessary to draw the system
        #inputs: list of link lengths, list of y values for each link at the given d, total height, number of links, given displacement
        #outputs: all the data needed to animate the system
        
        #get width of link
        W = .075*H
        
        #get x-spacing for links
        spx = [(i+1)*xspacing for i in range(n)]
        
        #get displacement angles of links
        theta = [np.arccos((L[j]-y[j])/L[j]) for j in range(n)]
        for j in range(n):
            if theta[j] > np.pi/2:
                theta[j] = np.pi/2
        
        #centers of endpoints of each link, bottom, middle, top
        centers = [[[spx[j],0],[spx[j]+(L[j]/2)*np.sin(theta[j]),(L[j]/2)*np.cos(theta[j])],[spx[j],L[j]-y[j]]] for j in range(n)]
        
        #angles for circles at ends of links
        #naming convention: tc = top circle, bc = bottom circle, tl = top link, bl = bottom link
        #top link
        phitctl = [np.linspace(theta[j],theta[j]+np.pi,20) for j in range(n)]
        phibctl = [np.linspace(theta[j]+np.pi,theta[j]+np.pi*2,20) for j in range(n)]
        #bottom link
        phitcbl = [np.linspace(-theta[j],-theta[j]+np.pi,20) for j in range(n)]
        phibcbl = [np.linspace(-theta[j]+np.pi,-theta[j]+2*np.pi,20) for j in range(n)]
        
        #make lists of points for link outlines
        #bottom link
        bcbl = [[] for j in range(n)]
        tcbl = [[] for j in range(n)]
        #top link
        bctl = [[] for j in range(n)]
        tctl = [[] for j in range(n)]
        
        for j in range(n):
            #bottom link
            bcbl[j] = [[(W/2)*np.cos(x)+centers[j][0][0],(W/2)*np.sin(x)+centers[j][0][1]] for x in phibcbl[j]]
            tcbl[j] = [[(W/2)*np.cos(x)+centers[j][1][0],(W/2)*np.sin(x)+centers[j][1][1]] for x in phitcbl[j]]
            #top link
            bctl[j] = [[(W/2)*np.cos(x)+centers[j][1][0],(W/2)*np.sin(x)+centers[j][1][1]] for x in phibctl[j]]
            tctl[j] = [[(W/2)*np.cos(x)+centers[j][2][0],(W/2)*np.sin(x)+centers[j][2][1]] for x in phitctl[j]]
            
        
        
        #make lists of vertices and codes for the path function to draw the links
        #start making vertices
        bot_verts = bcbl[0] + tcbl[0] + [bcbl[0][0]]
        top_verts = bctl[0] + tctl[0] + [bctl[0][0]]
        #start making codes
        codes = [Path.LINETO]*len(top_verts)
        codes[0] = Path.MOVETO
        codes[len(codes)-1] = Path.CLOSEPOLY
        
        #make the rest of the vertices
        for j in range(1,n):
            bot_verts = bot_verts + bcbl[j] + tcbl[j] + [bcbl[j][0]]
            top_verts = top_verts + bctl[j] + tctl[j] + [bctl[j][0]]
            
        #make the rest of the codes (repeat pattern for each link)
        codes = codes*n
        
        #make paths
        top_path = Path(top_verts,codes)
        bot_path = Path(bot_verts,codes)
        
        #make springs
        #initialize
        springx = [[] for j in range(n)]
        springy = [[] for j in range(n)]
        #calculate
        for j in range(n):
            #find number of spring section centers
            nc = int(np.ceil((H-L[j])/W))
            #find remainder to see if there's an extra bend or not
            rem = ((H-L[j])/W)%1
            #find length of leftover bit
            lb0 = H-L[j]-W*(nc-1)
            lbfrac = lb0/(H-L[j])
            lb = lbfrac*(H-L[j]-d+y[j])
            #find new distance between centers
            D = (H-L[j]-d+y[j]-lb)/(nc-1)
            #set spring bottom endpoints
            springx[j].append(centers[j][2][0])
            springy[j].append(centers[j][2][1])
            #find the rest of the bend locations
            for b in range(1,nc):
                sign_x = (-1)**(b-1)
                sx = centers[j][2][0] + sign_x*W*np.cos(np.arcsin(D/(2*W)))
                if b == 1:
                    sy = springy[j][b-1] + D/2
                else:
                    sy = springy[j][b-1] + D
                springx[j].append(sx)
                springy[j].append(sy)
            #find top endpoints
            sign_x = (-1)**(nc+1)
            #if no extra bend
            if rem<=.5:
                sx = centers[j][2][0] + sign_x*(lb/np.tan(np.arcsin(D/(2*W))))
                sy = H-d
                springx[j].append(sx)
                springy[j].append(sy)
            #if extra bend
            else:
                sx = centers[j][2][0] + sign_x*W*np.cos(np.arcsin(D/(2*W)))
                sy = springy[j][nc-1] + D
                springx[j].append(sx)
                springy[j].append(sy)
                sx = sx - sign_x*((lb-D/2)/np.tan(np.arcsin(D/(2*W))))
                sy = H-d
                springx[j].append(sx)
                springy[j].append(sy)    
        #add break between springs
        spr1x = []
        spr1y = []
        for j in range(n):
            spr1x = spr1x + springx[j] + [float('nan')]
            spr1y = spr1y + springy[j] + [float('nan')]    
            
        
        #get bar coordinates
        xybar = ([0,xspacing*(n+1)],[H-d,H-d])
        
        #make dots on link joints
        dotsx = []
        dotsy = []
        for j in range(n):
            dotsx = dotsx + [centers[j][i][0] for i in range(3)]
            dotsy = dotsy + [centers[j][i][1] for i in range(3)]
        
        dots = [dotsx]+[dotsy]  
        
        #draw dots for masses if dynamics is on
        if dyn:
            mx = [centers [j][2][0] for j in range(n)]
            my = [centers [j][2][1] for j in range(n)]
            mass = [mx]+[my]
        #leave masses blank if dynamics is off
        else:
            mx = []
            my = []    
            mass = [mx]+[my]
                
        
        return bot_path, top_path, xybar, spr1x, spr1y, dots, mass     
        
        
    
    #Set up the figure for the animation
    
    #get x and y axis limits for the link animation
    xmin = -.1*H
    xmax = xspacing*(n+1) +.1*H
    if 7*(xmax-xmin)/20 < H*1.2:
        xmin = ((n+1)*xspacing)/2-1.8*H
        xmax = xmin + 3.6*H
    ymin = -.1*H
    ymax = ymin + 7*(xmax-xmin)/20
    
    #initialization
    
    #initialize figure
    fig = plt.figure()
    
    #set axes
    ax1 = fig.add_subplot(211)
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)
    ax2 = fig.add_subplot(212)
    ax2.set_xlim(0,H+.1*H)
    #set P vs d plot y axis limits depending on model type
    #for dynamics
    if dyn:
        ax2.set_ylim(min(P_avg)*1.1,max(P_avg)*1.1)
    #for quasistatic
    else:
        ax2.set_ylim(0,F_crit_avg*n*1.1)
    
    #initialize patch vertices and move codes
    verts = [[-10.,-10.],[-10.,-9.],[-9.,-9.],[-9.,-10.],[-10.,-10.]]
    codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]
    
    #set path
    path = Path(verts,codes)
    
    #create blank link patches
    top_patch = patches.PathPatch(path,facecolor='white',edgecolor='white')
    bot_patch = patches.PathPatch(path,facecolor='white',edgecolor='white')
    
    #create screen clearer
    clr = patches.Rectangle((-H*.6,-.1*H),H*n+.6*H+H*.6,H*(n)/2+.1*H+.1*H,facecolor='white',edgecolor='white')
    
    
    #draw initial (blank) lines and patches
    lineb, = ax1.plot([],[],'k',lw=7)
    linesp, = ax1.plot([],[],'blue',lw=2)
    dots, = ax1.plot([],[],'k.')
    mass, = ax1.plot([],[],'ro')
    ax1.add_patch(top_patch)
    ax1.add_patch(bot_patch)
    linePd, = ax2.plot([],[])
    linePd_avg, = ax2.plot([],[],'r--')
    
    #label figure
    ax2.legend(["Actual Force Value","Force Value for Average Length"],fontsize='small')
    ax1.set_title('Link Spring System')
    ax2.set_title('Force vs. Displacement')
    ax2.set_xlabel('Displacement (m)')
    ax1.set_ylabel('m')
    ax2.set_ylabel('Force (N)')
    
        
    
    def init():
        #initialization function for matplotlib animation
        ax1.add_patch(clr)
        #empty lines
        linesp.set_data([],[])
        lineb.set_data([],[])
        dots.set_data([],[])
        mass.set_data([],[])
        linePd.set_data([],[])
        linePd_avg.set_data([],[])
        #blank patches
        ax1.add_patch(top_patch)
        ax1.add_patch(bot_patch)
        #return all objects to be animated
        return clr,lineb,bot_patch,top_patch,linesp,dots,mass,linePd_avg,linePd
    
    #animate
    def animate(i):
        #clear frame
        clr = patches.Rectangle((-H*.6,-.1*H),H*n+.6*H+H*.6,H*(n)/2+.1*H+.1*H,facecolor='white',edgecolor='white')
        ax1.add_patch(clr)
        #get frame data from draw function
        bot_path,top_path,xybar,spr1x,spr1y,dot,masses = draw(L,y[frameskip*i],H,n,d[frameskip*i],xspacing)
        #set top bar data
        lineb.set_data(xybar[0],xybar[1])
        #set spring data
        linesp.set_data(spr1x,spr1y)
        #set patches for links
        top_patch = patches.PathPatch(top_path,facecolor='gray',lw=2)
        bot_patch = patches.PathPatch(bot_path,facecolor='gray',lw=2)
        ax1.add_patch(bot_patch)
        ax1.add_patch(top_patch)
        #add dots
        dots.set_data(dot[0],dot[1])
        #add masses
        mass.set_data(masses[0],masses[1])
        #plot force vs. displacement
        linePd.set_data(d[:frameskip*i],P[:frameskip*i])
        #plot average force vs. displacement
        linePd_avg.set_data(d[:frameskip*i],P_avg[:frameskip*i])
        return clr,lineb,bot_patch,top_patch,linesp,dots,mass,linePd_avg,linePd
    
    #get interval between frames
    if dyn:
        inter = t_step*1000
        if save_ani:
            inter = 100
    else:
        inter = (1/fps)*1000
        if save_ani:
            inter = 100
    
    #animate
    anim = animation.FuncAnimation(fig,animate,init_func=init,frames=(len(d)-1)/frameskip,interval=inter,blit=True,repeat=False)
    plt.show()
    
    return anim   
        