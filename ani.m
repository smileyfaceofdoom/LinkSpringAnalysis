function ani(L,y,H,n,d)
%This function draws the link(s) at each position.
%input: Array of link lengths, array of corresponding y values, max height,
%number of links
%output: drawing of the links at the given y value

%find width of links
W = .05*H;

for j = 1:n
    %find angle from vertical
    theta = acos((L(j) - y(j))/L(j));
    %find endpoints of links
    C_1(j,1) = .5 + H*(j-1);
    C_1(j,2) = L(j) - y(j);
    C_2(j,1) = C_1(j,1) + (L(j)/2)*sin(theta);
    C_2(j,2) = (L(j)/2)*cos(theta);
    C_3(j,1) = C_1(j,1);
    C_3(j,2) = 0;
    
    %find corners of the top link
    %x coordinates
    cor_tx(j,1) = C_1(j,1) + (W/2)*cos(theta);
    cor_tx(j,2) = C_1(j,1) - (W/2)*cos(theta);
    cor_tx(j,3) = C_2(j,1) - (W/2)*cos(theta);
    cor_tx(j,4) = C_2(j,1) + (W/2)*cos(theta);
    %y coordinates
    cor_ty(j,1) = C_1(j,2) + (W/2)*sin(theta);
    cor_ty(j,2) = C_1(j,2) - (W/2)*sin(theta);
    cor_ty(j,3) = C_2(j,2) - (W/2)*sin(theta);
    cor_ty(j,4) = C_2(j,2) + (W/2)*sin(theta);
    
    %find corners of the bottom link
    %x coordinates
    cor_bx(j,1) = C_2(j,1) - (W/2)*cos(theta);
    cor_bx(j,2) = C_2(j,1) + (W/2)*cos(theta);
    cor_bx(j,3) = C_3(j,1) + (W/2)*cos(theta);
    cor_bx(j,4) = C_3(j,1) - (W/2)*cos(theta);
    %y coordinates
    cor_by(j,1) = C_2(j,2) + (W/2)*sin(theta);
    cor_by(j,2) = C_2(j,2) - (W/2)*sin(theta);
    cor_by(j,3) = C_3(j,2) - (W/2)*sin(theta);
    cor_by(j,4) = C_3(j,2) + (W/2)*sin(theta);
    
    %create circles
    t = 2*(pi/1000)*(1:1000);
    r = (W/2);
    for i = 1:1000
    circ_1x(j,i) = C_1(j,1) + r*cos(t(i));
    circ_1y(j,i) = C_1(j,2) + r*sin(t(i));
    circ_2x(j,i) = C_2(j,1) + r*cos(t(i));
    circ_2y(j,i) = C_2(j,2) + r*sin(t(i));
    circ_3x(j,i) = C_3(j,1) + r*cos(t(i));
    circ_3y(j,i) = C_3(j,2) + r*sin(t(i));
    end
    %create spring
    %number of centers
    nc = ceil((H-L(j))/W);
    %set spring endpoint
    s_x(j,1) = C_1(j,1);
    s_y(j,1) = C_1(j,2);
    D = (W-(d-y(j))/(nc));
    %find locations of the bends
    for b = 2:nc+1
        sign_x = (-1)^b;
        s_x(j,b) = C_1(j,1) + sign_x*W*cos(asin(D/(2*W)));
        s_y(j,b) = C_1(j,2) + D*(b-1) - D/2;
    end
    lb = (H-d) - s_y(j,nc+1);
    s_x(j,nc+2) = C_1(j,1) + sign_x*(W*cos(asin(D/(2*W))) - lb/tan(asin(D/(2*W))));
    s_y(j,nc+2) = H-d;
    
end
%draw top block
block_x(1) = .25;
block_x(2) = n*H - .25;
block_x(3) = n*H - .25;
block_x(4) = .25;
block_y(1) = H-d;
block_y(2) = H-d;
block_y(3) = H-d + W/2;
block_y(4) = H-d + W/2;

%draw rectangles and circles
fill(block_x,block_y,'k')
hold on

for j = 1:n
    
    fill(cor_tx(j,:),cor_ty(j,:),'k')
    fill(cor_bx(j,:),cor_by(j,:),'k')
end
for j = 1:n
    fill(circ_1x(j,:),circ_1y(j,:),'k')
    fill(circ_2x(j,:),circ_2y(j,:),'k')
    fill(circ_3x(j,:),circ_3y(j,:),'k')
    
end
%draw springs
for j = 1:n
    plot(s_x(j,:),s_y(j,:),'k')
end
axis([0,(n*H+.5), (-.1*H), (H+.1)])
axis('equal')
hold off


end