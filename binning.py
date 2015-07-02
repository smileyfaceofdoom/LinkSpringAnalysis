def binning(r_p,i_max):
    import math
    #This function applies binning to a set of random numbers to reduce computational cost. It uses 2 Gauss points per cell
    #inputs: the list of points to be binned (r_p), the number of bins (i_max)
    #outputs: The binned points (r_b), the weights of the binned points (W_b), the new number of points (n)
    
    #set number of points and nodes 
    p_max = len(r_p) #number of points
    
    #find node locations, create r_i
    r_i1 = range(i_max)
    r_i = [x/float(i_max-1) for x in r_i1]
    
    
    #assign points to the cells they're in
    num_cells = i_max-1
    #initialize cell list
    cells = []
    #for each cell
    for j in range(num_cells):
        cells.append([]) #create cell
        #for each data point
        for k in range(p_max):
            #check to see if data point is in cell
            if r_p[k] >= r_i[j] and r_p[k] <r_i[j+1]:
                #if so, assign to cell
                cells[j].append(r_p[k])
    
    
    #initialize nodal volume and point count to zero for all nodes
    V = [0]*i_max
    N = [0]*i_max
    
    #initialize counter
    f = 1
    #for each cell:
    while f<i_max:
        #find nodes on either side of the cell
        a = r_i[f-1]
        b = r_i[f]
        #find cell length
        Lc = b-a
        
        #find cell's contributions to nodal volume
        V[f-1] += Lc/2
        V[f] += Lc/2
        
        #find cell's contributions to nodal point count for each point in cell
        for j in range(len(cells[f-1])):
            #evaluate shape functions at the point r_p[j]
            S_1 = (cells[f-1][j]-b)/(a-b)
            S_2 = (cells[f-1][j]-a)/(b-a)
            
            #add to nodal point count
            N[f-1] += S_1
            N[f] += S_2
        
        f+=1 #increment counter 
    
    #calculate nodal weights
    w = [Ni/p_max for Ni in N]
    
    #calculate occupancies
    f=1 #initialize counter variable
    r_g = [] #initialize gauss point locations
    W_g = [] #initialize occupancies
    #for each cell
    while f<i_max:
        #find nodes on either side of cell
        a = r_i[f-1]
        b = r_i[f]
        Lc = b-a
        
        #set gauss point locations
        gp1 = (r_i[f-1] + r_i[f])/2 - (Lc/2)*(1/math.sqrt(3))
        gp2 = (r_i[f-1] + r_i[f])/2 + (Lc/2)*(1/math.sqrt(3))
        r_g.append(gp1)
        r_g.append(gp2)
        
        #set standard Gauss weights
        gw_1 = 1
        gw_2 = 1
        
        #find scaled Gauss weights
        k_g1 = gw_1*(b-a)/2
        k_g2 = gw_2*(b-a)/2
        
        #evaluate shape functions at Gauss points
        S1_g1 = (gp1-b)/(a-b)
        S2_g1 = (gp1-a)/(b-a)
        S1_g2 = (gp2-b)/(a-b)
        S2_g2 = (gp2-a)/(b-a)
        
        #find delta functions
        delta1_rg1 = S1_g1/V[f-1]
        delta2_rg1 = S2_g1/V[f]
        delta1_rg2 = S1_g2/V[f-1]
        delta2_rg2 = S2_g2/V[f]
        
        #find q(r_g)
        q_rg1 = w[f-1]*delta1_rg1 + w[f]*delta2_rg1
        q_rg2 = w[f-1]*delta1_rg2 + w[f]*delta2_rg2
        
        #find occupancies
        W_g1 = k_g1*q_rg1
        W_g2 = k_g2*q_rg2
        
        W_g.append(W_g1)
        W_g.append(W_g2)
        
        
        f+=1 #increment counter  
        
    #remove occupancies of zero
    W_b = []
    r_b = []
    for Wg, rg in zip(W_g,r_g):
        if Wg != 0:
            W_b.append(Wg)
            r_b.append(rg)
    
    #find new number of data points
    n = len(r_b)
            
    return W_b,r_b,n                       
    

