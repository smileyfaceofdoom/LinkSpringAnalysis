%Katharin Jensen
%This code creates a force vs. displacement plot for a model [describe
%model better later]
tic

clear
clc
close all

%initialize
ea_t = 10; %other spring constant type thing, N
k_t = 100; %spring constant, N/m
H = 1; %max height, m
%n = 100; %number of links
m = 4;
n_l = 10;
i_max = 11;
ea = ea_t/n_l;
k = k_t/n_l;

%create array of lengths
%generate random numbers
%R = rand(1,n_l);
R = linspace(.3,.9,n_l);
%[W_b,r_b,n] = binning(R,i_max);
n=n_l;
%create Weibull distribution
l = H.*((log(1./R)/log(2)).^(1/m));
%create exponential distribution of lengths
L = H.*(1-exp(-l));
%L = H*[.7,.3];
%find average length
l_avg = H.*((log(1./R)/log(2)).^(1/m));
Lexp_avg = H.*(1-exp(-l_avg));
L_avg = mean(Lexp_avg);

%create displacement distribution
d = linspace(0,H,1000);
%indicator for whether each link has buckled. 0 = not yet, 1 = yes, 2 =
%completely collapsed
c = zeros(n,1);
c_avg = 0;
F_crit = k.*L/4; %calculate array of critical buckling values
F_old = F_crit; %initialize fzero guesses
F_crit_avg = k.*L_avg/4; %calculate critical buckling value for average
display(F_crit_avg*n_l)
F_old_avg = F_crit_avg; %initialize fzero guess for average
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
%     P(i) = 0; %initialize force
%     for j = 1:n
%         P(i) = P(i) + F(j)*W_b(j)*n_l;
%     end
    
    %find force for average length
    if (c_avg == 0 && d(i) >= (H - L_avg))
            c_avg = 1;
    end
    %if link hasn't buckled
    if (c_avg == 0)
        %calculate link force
        F_avg = ea*atanh(d(i)/(H-L_avg));
        %check if link has buckled
        if (F_avg >= F_crit_avg)
            c_avg = 1; %set indicator if it has
        end
    %if link has buckled but not completely collapsed    
    elseif (c_avg == 1)
        %numerically find the link force
        F_avg = fzero(@(f) L_avg*(1 - (4*f)/(k*L_avg)) + (H-L_avg)*tanh(f/ea) - d(i),F_old_avg);
        F_old_avg = F_avg;
        %check if link has completely collapsed
        %compute y
        y_avg = L_avg - (4*F_avg)/k;
        if (y_avg >= L_avg)
            c_avg = 2; %set indicator if it has
            F_avg = ea*atanh((d(i)-L_avg)/(H-L_avg));
        end
    %once the link has completely collapsed    
    else
        %calculate link force
        F_avg = ea*atanh((d(i)-L_avg)/(H-L_avg));
    end    
    P_avg(i) = n_l*F_avg;
end

display(L)
%axis([0,1,-10,200])

%IMPORTANT NOTE: The below animations CANNOT be used if binning was as
%well, because it requires the link displacement values for each individual
%link, which you do NOT get if binning is applied.

%---Uncomment the below section to make a nice animation of the system.
%Note, it's really slow, so don't use it for more than a few links. The
%double-commented code is specifically for making a .gif of the animation,
%so only completely uncomment it if that's what you want to do.---

%figure
% %fn = 'TestAni2.gif'; %filename for making a .gif, change to whatever
% for i = 1:1000
%     subplot(2,1,1);
%     %figure(2)
%     ani(L,y(:,i),H,n,d(i));
%     subplot(2,1,2);
%     %figure(3)
%     plot(d(1:i),P(1:i));
%     title('Force vs. Displacement')
%     ylabel('Force (N)')
%     xlabel('Displacement (m)')
%     axis([0,H+.1,0,10*n+15]) %y axis scale should probably be adjusted
%     drawnow
%
% %    %make .gif of animation
% %    M = getframe(gcf);
% %    im = frame2im(M);
% %    [A,map] = rgb2ind(im,256);
% %    if i==1
% %        imwrite(A,map,fn,'gif','LoopCount',1,'DelayTime',.05);
% %    else
% %        imwrite(A,map,fn,'gif','WriteMode','append','DelayTime',.05);
% %    end
% 
% end

%---End nice animation section---

%---Uncomment the code below to make a rough animation of the system. It's
%faster than the nice animation, so it can handle more links, but it gets
%slow (and tiny) after about 10, so try not to do many more than that. The
%double-commented code makes a .gif of the animation, so don't completely
%uncomment that part unless that's what you want to do.---

% %fn = 'Rough10Link.gif'; %filename for making a .gif, change to whatever
% for i = 1:1000
%     subplot(2,1,1);
%     %figure(2)
%     ani_rough(L,y(:,i),H,n,d(i));
%     subplot(2,1,2);
%     %figure(3)
%     plot(d(1:i),P(1:i));
%     title('Force vs. Displacement')
%     ylabel('Force (N)')
%     xlabel('Displacement (m)')
%     axis([0,H+.1,0,F_crit_avg*n+20])
%     drawnow
%     
% %    %make .gif of animation
% %    M = getframe(gcf);
% %    im = frame2im(M);
% %    [A,map] = rgb2ind(im,256);
% %    if i==1
% %        imwrite(A,map,fn,'gif','LoopCount',1,'DelayTime',.05);
% %    else
% %        imwrite(A,map,fn,'gif','WriteMode','append','DelayTime',.05);
% %    end
% 
% end

%---End rough animation section---

figure
plot(d,P)
hold on
plot(d,P_avg,'r--');
hold off
title('Force vs. Displacement')
xlabel('Displacement (m)')
ylabel('Force (N)')
legend('Actual Value','Value With Average Length')
axis([0, H+.1*H, 0, F_crit_avg*n_l+20])

toc






