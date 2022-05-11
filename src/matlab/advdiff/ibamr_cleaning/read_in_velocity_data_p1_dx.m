function read_in_velocity_data_p1()
%function read_in_velocity_data_p1.m 
%input: 
%output: 
%sets xlength and ylength based on the experimental piv data

%this function sets it so that the computational 
%domain for c sits inside the piv exp data

global xlength ylength Nx Ny dx dy
global piv_data_filename piv_data_filename_interior
global xshift_piv_data yshift_piv_data
global coarsedx coarsedy 


%filename
load(['pivdata/' piv_data_filename '.mat'])   

%setting x_length and y_length (data is in cm -> need mm)
xpiv = eval([piv_data_filename_interior.filename '.' piv_data_filename_interior.x]);      
ypiv = eval([piv_data_filename_interior.filename '.' piv_data_filename_interior.y]);

xlength = 10*(min(xpiv(:,end))-max(xpiv(:,1)))
ylength = 10*(min(ypiv(end,:))-max(ypiv(1,:)))

Nxcoarse = floor(xlength/(coarsedx));
Nycoarse = floor(ylength/(coarsedy));

%minus 4 allows for extrap for weno2
xlength = (Nxcoarse-4)*coarsedx
ylength = (Nycoarse-4)*coarsedy

xshift_piv_data = -10*max(xpiv(:,1))-2*coarsedx;
yshift_piv_data = -10*max(ypiv(1,:))-2*coarsedy;

%for convergence tests this is better to set manually or as set above
% xlength = 2.1250;
% ylength = 1.8750;
% 
% Nx = round(xlength/dx);
% Ny = round(ylength/dy);

pause











