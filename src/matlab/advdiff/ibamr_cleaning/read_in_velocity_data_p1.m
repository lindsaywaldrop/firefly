function [parameters, simulation] = read_in_velocity_data_p1(paths, parameters, simulation)
%function read_in_velocity_data_p1.m 
%input: 
%output: 
%sets xlength and ylength based on the experimental piv data

%this function sets it so that the computational 
%domain for c sits inside the piv exp data

%NOTES: SK 2017_11_16
%COMMENTED ALL RETURN STUFF OUT 

%global xlength ylength Nx Ny dx dy pathbase_piv
%global piv_data_filename piv_data_filename_interior
%global xshift_piv_data yshift_piv_data
%global Nxcoarse 

strcat(paths.pathbase_ibamr, parameters.ibamr_data_filename, '.mat')

%filenames
flickdata = load([paths.pathbase_ibamr parameters.ibamr_data_filename '.mat']);

%finding the endpts for x_length and y_length (data is in m) 
flick_xpiv = eval(['flickdata.' parameters.ibamr_data_filename_interior.x]);     
flick_ypiv = eval(['flickdata.' parameters.ibamr_data_filename_interior.y]);

parameters.domainlimits = set_domain_limits(paths, parameters, simulation.run_id);

xLmin = parameters.domainlimits(1); 
xLmax = parameters.domainlimits(2); 
yLmin = parameters.domainlimits(3); 
yLmax = parameters.domainlimits(4); 

%setting x_length and y_length 

xlength_newd = xLmax - xLmin; 
ylength_newd = yLmax - yLmin; 


%coarsedx = (max(max(flick_xpiv))-min(min(flick_xpiv)))/parameters.Nxcoarse;
%coarsedy = coarsedx;

%parameters.Nycoarse = floor(ylength/coarsedy);
%parameters.Nycoarse_shift = coarsedy*(ylength/coarsedy - floor(ylength/coarsedy))/2;

parameters.dx = (max(max(flick_xpiv))-min(min(flick_xpiv)))/parameters.GridSize;
parameters.dy = (max(max(flick_ypiv))-min(min(flick_ypiv)))/parameters.GridSize;

parameters.Nx = floor(xlength_newd/parameters.dx);
parameters.Ny = floor(ylength_newd/parameters.dy);

%minus 4 allows for extrap for weno2
parameters.xlength = (parameters.Nx-4)*parameters.dx;
parameters.ylength = (parameters.Ny-4)*parameters.dy;

%flickdata shift 
%parameters.xshift_ibamr_data(1) = -xLmin-2*coarsedx;
%parameters.yshift_ibamr_data(1) = -yLmin-2*coarsedy-parameters.Nycoarse_shift;

parameters.xshift_ibamr_data(1) = 0;
parameters.yshift_ibamr_data(1) = 0;


