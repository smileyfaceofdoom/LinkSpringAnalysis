%Katharin Jensen
%This code creates a force vs. displacement plot for a model [describe
%model better later]
tic

clear
clc
close all

%initialize
ea = 10; %other spring constant type thing, N
k = 100; %spring constant, N/m
H = 1; %max height, m
n = 2; %number of links
m = 4;

%create array of lengths
%generate random numbers
R = rand(1,n);
%R = [.1,.2,.3,.4,.5,.6,.7,.8,.9];

%create Weibull distribution
l = H.*((log(1./R)/log(2)).^(1/m));
%verification


%create exponential distribution of lengths
L = H.*(1-exp(-l));
%L = H*[.3,.5001];

%create displacement distribution
d = linspace(0,H,1000);
%indicator for whether each link has buckled. 0 = not yet, 1 = yes, 2 =
%completely collapsed
c = zeros(n,1);
F_crit = k.*L/4; %calculate array of critical buckling values
F_old = F_crit; %initialize fzero guesses
y = 0; %initial link displacement
%solve for P
for i = 1:1000
    %count iterations
    %z = i;
    %for each link
    for j = 1:n
        %check if displacement has reached the link
        if (c(j) == 0 && d(i) >= (H - L(j)))
            c(j) = 1;
        end
        %if link hasn't buckled
        if (c(j) == 0)
            %calculate link force
            F(j) = ea*atanh(d(i)/(H-L(j)));
            y(j,i) = 0;
            %check if link has buckled
            if (F(j) >= F_crit(j))
                c(j) = 1; %set indicator if it has
                crit_it(j) = i;
            end
        %if link has buckled but not completely collapsed    
        elseif (c(j) == 1)
            %numerically find the link force
            F(j) = fzero(@(f) L(j)*(1 - (4*f)/(k*L(j))) + (H-L(j))*tanh(f/ea) - d(i),F_old(j));
            F_old(j) = F(j);
            %check if link has completely collapsed
            %compute y
            y(j,i) = L(j) - (4*F(j))/k;
            if (y(j,i) >= L(j))
                c(j) = 2; %set indicator if it has
                F(j) = ea*atanh((d(i)-L(j))/(H-L(j)));
            end
        %once the link has completely collapsed    
        else
            %calculate link force
            F(j) = ea*atanh((d(i)-L(j))/(H-L(j)));
            y(j,i) = L(j);
        end
    end

    %sum all forces to get the total force    
    P(i) = sum(F);
end

display(L)
plot(d,P)
title('Force vs. Displacement')
xlabel('Displacement (m)')
ylabel('Force (N)')
%axis([0,1,-10,200])

%figure
%make animation
for i = 1:1000
    subplot(2,1,1);
    %figure(2)
    ani(L,y(:,i),H,n,d(i));
    subplot(2,1,2);
    %figure(3)
    plot(d(1:i),P(1:i));
    title('Force vs. Displacement')
    ylabel('Force (N)')
    xlabel('Displacement (m)')
    axis([0,H+.1,0,10*n+15]) %y axis scale should probably be adjusted
    drawnow
end


toc




