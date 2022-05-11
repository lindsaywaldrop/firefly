load('viz_IB2d1.mat')  % Loads interpolated and non-interpolated data from VisIt
load('hairinfo.mat')   % Loads point boundary point data

%w=magnitude(x,y,Vinterp.U_x,Vinterp.U_y); % Calculate magnitudes from interpolated velocity data
%pcolor(x,y,w),shading flat  %pcolor plot
skip = 10; 
quiver(x(1:skip:end,1:skip:end),y(1:skip:end,1:skip:end),Vinterp.U_x(1:skip:end,1:skip:end),Vinterp.U_y(1:skip:end,1:skip:end));
hold on
plot(p.h1_x,p.h1_y,'k.')  %plot hair 1 (center hair)
plot(p.h2_x,p.h2_y,'k.')  %plot hair 2 (top hair)
plot(p.h3_x,p.h3_y,'k.')  %plot hair 3 (bottom hair)
hold off

%p.hdia - diameter of each hair 
%hair1Centerx - center of hair 1 - x position 
%hair1Centery - center of hair 1 - y position 


%meters and seconds 
%only flick 
%run to time 0.1 seconds 
% domain 4 by 2 amd resolution is not going to  change 
% initial conditions - stripe concentration 
%                      - continuous incoming concentration 
%start it close to the hairs and make it an adjustable parameter 

%Run convergence of this data set 
%1233 runs 
%aim 1-2 hours per run 

%domain size - think about making the domain smaller....

%Questions
%(1) Is it always going to be three hairs? Yes, mostly
%(2) Automating max velocity for dtfactor 
%(3) Automating domain size based on when velocity is modified from
%horizontal 
