"""This will animate the link spring system"""
def ani(L,y,H,n,d):
    
    #import
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import matplotlib.animation as animation
    import numpy as np
    
    #get width of link
    W = .05*H
    
    #get x-spacing for links
    spx = [H*i+H/2 for i in range(n)]
    
    #get displacement angles of links
    theta = [np.arccos((L[j]-y[j])/L[j]) for j in range(n)]
    for j in range(n):
        if theta[j] > np.pi/2:
            theta[j] = np.pi/2
    
    #centers of endpoints of each link, bottom, middle, top
    centers = [[[spx[j],0],[spx[j]+(L[j]/2)*np.sin(theta[j]),(L[j]/2)*np.cos(theta[j])],[spx[j],L[j]-y[j]]] for j in range(n)]
    
    #angles for circles at ends of links
    phitctl = [np.linspace(theta[j],theta[j]+np.pi,20) for j in range(n)]
    phibctl = [np.linspace(theta[j]+np.pi,theta[j]+np.pi*2,20) for j in range(n)]
    phitcbl = [np.linspace(-theta[j],-theta[j]+np.pi,20) for j in range(n)]
    phibcbl = [np.linspace(-theta[j]+np.pi,-theta[j]+2*np.pi,20) for j in range(n)]
    
    #make lists of points for link outlines
    bcbl = [[] for j in range(n)]
    tcbl = [[] for j in range(n)]
    bctl = [[] for j in range(n)]
    tctl = [[] for j in range(n)]
    
    for j in range(n):
        bcbl[j] = [[(W/2)*np.cos(x)+centers[j][0][0],(W/2)*np.sin(x)+centers[j][0][1]] for x in phibcbl[j]]
        tcbl[j] = [[(W/2)*np.cos(x)+centers[j][1][0],(W/2)*np.sin(x)+centers[j][1][1]] for x in phitcbl[j]]
        bctl[j] = [[(W/2)*np.cos(x)+centers[j][1][0],(W/2)*np.sin(x)+centers[j][1][1]] for x in phibctl[j]]
        tctl[j] = [[(W/2)*np.cos(x)+centers[j][2][0],(W/2)*np.sin(x)+centers[j][2][1]] for x in phitctl[j]]
        
    
    
    #make lists of verts and codes for the path function to draw the links
    #start making verts
    verts = bcbl[0] + tcbl[0] + [bcbl[0][0]]
    
    #start making codes
    codes = [Path.LINETO]*len(verts)
    codes[0] = Path.MOVETO
    codes[len(codes)-1] = Path.CLOSEPOLY
    
    #make the rest of the verts
    verts = verts + bctl[0] + tctl[0] + [bctl[0][0]]
    for j in range(1,n):
        verts = verts + bcbl[j] + tcbl[j] + [bcbl[j][0]] + bctl[j] + tctl[j] + [bctl[j][0]]
        
    #make the rest of the codes
    codes = codes*2*n
    
    path = Path(verts,codes)
    
    #make springs
    #initialize
    springx = [[] for j in range(n)]
    springy = [[] for j in range(n)]
    #calculate
    for j in range(n):
        #find number of spring section centers
        nc = int(np.ceil((H-L[j])/W))
        #find leftover bit
        rem = ((H-L[j])/W)%1
        lb0 = H-L[j]-W*(nc-1)
        lbfrac = lb0/(H-L[j])
        lb = lbfrac*(H-L[j]-d+y[j])
        #find amount each spring bend displaces
        D = (H-L[j]-d+y[j]-lb)/nc
        #set spring endpoints
        springx[j].append(centers[j][2][0])
        springy[j].append(centers[j][2][1])
        #find the rest of the bend locations
        for b in range(1,nc+1):
            sign_x = (-1)**(b-1)
            sx = centers[j][2][0] + sign_x*W*np.cos(np.arcsin(D/(2*W)))
            sy = springy[j][b-1] + D
            springx[j].append(sx)
            springy[j].append(sy)
        
        sign_x = (-1)**(nc+2)
        if rem<=.5:
            sx = centers[j][2][0] + sign_x*(lb/np.tan(np.arcsin(D/(2*W))))
            sy = H-d
            springx[j].append(sx)
            springy[j].append(sy)
        else:
            sx = centers[j][2][0] + sign_x*W*np.cos(np.arcsin(D/(2*W)))
            sy = springy[j][nc] + D
            springx[j].append(sx)
            springy[j].append(sy)
            sx = sx - sign_x*((lb-D/2)/np.tan(np.arcsin(D/(2*W))))
            sy = H-d
            springx[j].append(sx)
            springy[j].append(sy)    
    #add break in springs
    spr1x = []
    spr1y = []
    for j in range(n):
        spr1x = spr1x + springx[j] + [float('nan')]
        spr1y = spr1y + springy[j] + [float('nan')]    
        
    
    #get bar coordinates
    xybar = ([0,H*n],[H-d,H-d])
    
    #make dots on link joints
    dotsx = []
    dotsy = []
    for j in range(n):
        dotsx = dotsx + [centers[j][i][0] for i in range(3)]
        dotsy = dotsy + [centers[j][i][1] for i in range(3)]
    
    dots = [dotsx]+[dotsy]   
    
    
    '''fig = plt.figure()
    ax = fig.add_subplot(111)
    #draw links
    patch = patches.PathPatch(path,facecolor="gray",lw=2)
    ax.add_patch(patch)
    #draw top rectangle
    ax.plot([0,H*n],[H-d,H-d],'brown',lw=10)
    #draw springs
    ax.plot(spr1x,spr1y,'k',lw=2)
    ax.plot(dots[0],dots[1],'k.')
    ax.axis('equal')
    plt.show()'''     
    
    return path, xybar, spr1x, spr1y, dots     
    
    
    
    
    
    
    
    
    
ani([.62,.41],[.25,.3],1.0,2,.3)    
    