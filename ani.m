function ani(L,y,H,n,d)
%This function draws the link(s) at each position.
%input: Array of link lengths, array of corresponding y values, max height,
%number of links
%output: drawing of the links at the given y value

%find width of links
W = .05*H;

%initialize spring matrices
s_x = zeros(n,20);
s_y = zeros(n,20);

for j = 1:n
    %find angle from vertical
    theta = acos((L(j) - y(j))/L(j));
    if theta > pi/2
        theta = pi/2;
    end
    %find endpoints of links
    C_1(j,1) = .5*H + H*(j-1);
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
    %get spring coordinates
    
    %create spring
    %number of centers
    nc = ceil((H-L(j))/W);
    %indicator of where the top of the spring falls on the board
    rem = (H-L(j))/W - floor((H-L(j))/W);
    %find leftover if spring length doesn't divide evenly by width
    lb_init = (H-L(j)) - W*(nc-1);
    lb_frac = lb_init/(H-L(j));
    lb = lb_frac*((H-L(j))-d+y(j));
    %find amount each bend of spring displaces
    D = ((H-L(j))-d+y(j)-lb)/nc;
    %set spring endpoint
    s_x(j,1) = C_1(j,1);
    s_y(j,1) = C_1(j,2);
    %find locations of the bends
    for b = 2:nc+1
        sign_x = (-1)^b;
        s_x(j,b) = C_1(j,1) + sign_x*W*cos(asin(D/(2*W)));
        s_y(j,b) = C_1(j,2) + D*(b-1) - D/2;
    end
    sign_x = (-1)^(nc+2);
    if rem<.5
        s_x(j,nc+2) = C_1(j,1) + sign_x*(lb/tan(asin(D/(2*W))));
        s_y(j,nc+2) = H-d;
    else
        s_x(j,nc+2) = C_1(j,1) + sign_x*W*cos(asin(D/(2*W)));
        s_y(j,nc+2) = C_1(j,2) + D*(nc+1) - D/2;
        s_x(j,nc+3) = s_x(j,nc+2) - sign_x*((lb-D/2)/tan(asin(D/(2*W))));
        s_y(j,nc+3) = H-d;
    end

end
%draw top block
block_x(1) = .25*H;
block_x(2) = n*H - .25*H;
block_x(3) = n*H - .25*H;
block_x(4) = .25*H;
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
    s_x_j = s_x(j,:);
    %s_y_j = s_y(j,:);
    ind = find(s_x_j == 0,1);
    %indy = find(s_y_j ==0,1);
    %disp(indx)
    %disp(indy)
    %disp(d)
    %disp(y)
    plot(s_x(j,1:(ind-1)),s_y(j,1:(ind-1)))
end
axis([0,(n*H+.5*H), (-.1*H), (H+.1*H)])
axis('equal')
hold off


end