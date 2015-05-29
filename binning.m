function [W_b,r_b,n] = binning(r_p,i_max)
%This function applies binning to the set of random numbers in order to
%reduce the computational cost. It uses 2 Gauss points per cell.

%set number of points and nodes
p_max = length(r_p); %number of points

%find node locations, create r_i
r_i=zeros(1,i_max);
for j = 2:i_max
    r_i(j) = r_i(j-1) + 1/(i_max-1);
end

%assign points to the cells they're in
%initialize matrix, made of values r_p can't have. Rows correspond to cells
%and columns to data points.
cells = ones((i_max-1),p_max);
c = ones(size(cells,1),1); %initialize column counter for each cell
%for each cell
for j = 1:size(cells,1)
    
    %for each data point
    for k = 1:p_max
        
        %if the data point is in the cell
        if r_p(k)>=r_i(j) && r_p(k)<r_i(j+1)
            
            cells(j,c(j)) = r_p(k); %assign point to cell
            c(j) = c(j) + 1; %increment column counter to next column
        
        end
        
    end
    
end

%initialize nodal volume and point count both to zero for all nodes
V = zeros(1,i_max);
N = zeros(1,i_max);

f = 1; %initialize counter variable
%for each cell
while f<i_max
    
    %find nodes on either side of the cell
    a = r_i(f);
    b = r_i(f+1);
    %find cell length
    L = b-a;
    
    %find cell's contributions to nodal volume
    V(f) = V(f) + L/2;
    V(f+1) = V(f+1) + L/2;
    
    %find cell's contributions to nodal point count
    %for each point in the cell (i.e. only the columns that had a point
    %value put in them, noting that the value of c(f) is one greater than
    %that number of columns)
    for j = 1:(c(f)-1)
        
        %evaluate shape functions at the point r_p(j)
        S_1 = (cells(f,j)-b)/(a-b);
        S_2 = (cells(f,j)-a)/(b-a);
        
        %add to nodal point count
        N(f) = N(f) + S_1;
        N(f+1) = N(f+1) + S_2;
        
    end
    
    f = f+1;
    
end

%calculate nodal weights
w = N./p_max;

%calculate occupancies
f = 1; %re-initialize counter variable
%for each cell
while f<i_max
    
    %find nodes on either side of the cell
    a = r_i(f);
    b = r_i(f+1);
    L = b-a;
    
    %set Gauss point locations
    r_g(2*f-1) = (r_i(f)+r_i(f+1))/2 - (L/2)*(1/sqrt(3));
    r_g(2*f) = (r_i(f)+r_i(f+1))/2 + (L/2)*(1/sqrt(3));
    
    %set standard Gauss weights
    gw_1 = 1;
    gw_2 = 1;
    
    %find scaled Gauss weights
    k_g1 = gw_1*(b-a)/2;
    k_g2 = gw_2*(b-a)/2;
    
    %evaluate shape functions at Gauss points
    S1_g1 = (r_g(2*f-1)-b)/(a-b);
    S2_g1 = (r_g(2*f-1)-a)/(b-a);
    S1_g2 = (r_g(2*f)-b)/(a-b);
    S2_g2 = (r_g(2*f)-a)/(b-a);
    
    %find delta functions
    delta1_rg1 = S1_g1/V(f);
    delta2_rg1 = S2_g1/V(f+1);
    delta1_rg2 = S1_g2/V(f);
    delta2_rg2 = S2_g2/V(f+1);
    
    %find q(r_g)
    q_rg1 = w(f)*delta1_rg1 + w(f+1)*delta2_rg1;
    q_rg2 = w(f)*delta1_rg2 + w(f+1)*delta2_rg2;

    
    %find occupancies
    W_g(2*f-1) = k_g1*q_rg1;
    W_g(2*f) = k_g2*q_rg2;
    
    f = f+1;
    
end

%discard all occupancies of zero
z = 1; %initialize counter to make sure X_b and W_b don't have skipped indeces
for j = 1:length(W_g)
    
    if W_g(j) ~= 0
        
        W_b(z) = W_g(j);
        r_b(z) = r_g(j);
        z = z+1;
        
    end
    
end
n = length(r_b);

end