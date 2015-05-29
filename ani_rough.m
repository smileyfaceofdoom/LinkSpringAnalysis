function ani_rough(L,y,H,n,d)
%This function generates a rough image of the links and springs as just
%lines, for the sake of being able to quickly see the deformation of a
%large number of links.

%create coordinate matrices for links and springs
%initialize matrices
coord_x = zeros(n,4);
coord_y = zeros(n,4);
for j = 1:n
    
    %find angle from vertical
    theta = acos((L(j) - y(j))/L(j));
    if theta > pi/2
        theta = pi/2;
    end
    
    %find ends of links
    coord_x(j,1) = .5*H + H*(j-1);
    coord_y(j,1) = 0;
    coord_x(j,2) = coord_x(j,1) + (L(j)/2)*sin(theta);
    coord_y(j,2) = (L(j)/2)*cos(theta);
    coord_x(j,3) = coord_x(j,1);
    coord_y(j,3) = L(j) - y(j);
    
    %find top end of spring
    coord_x(j,4) = coord_x(j,3);
    coord_y(j,4) = H - d;
    
end

%find ends of top bar
bar(1,1) = 0;
bar(1,2) = H*n;
bar(2,1) = H-d;
bar(2,2) = H-d;

%plot bar
plot(bar(1,:),bar(2,:),'k','LineWidth',4)
hold on
%plot links and springs
for j = 1:n
    %plot links
    plot(coord_x(j,1:3),coord_y(j,1:3),'k','LineWidth',2)
    %plot springs
    plot(coord_x(j,3:4),coord_y(j,3:4),'m')
end
axis([0,(n*H+.5*H), (-.1*H), (H+.1*H)])
axis('equal')
hold off


end