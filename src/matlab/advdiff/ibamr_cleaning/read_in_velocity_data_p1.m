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

strcat(paths.pathbase_piv, parameters.piv_data_filename, '.mat')

%filenames
flickdata = load([paths.pathbase_piv parameters.piv_data_filename '.mat']);

%finding the endpts for x_length and y_length (data is in m) 
flick_xpiv = eval(['flickdata.' parameters.piv_data_filename_interior.x]);     
flick_ypiv = eval(['flickdata.' parameters.piv_data_filename_interior.y]);


if  strcmp(parameters.domainlimits,'auto')
    xLmin = parameters.piv_data_filename_interior.conversion_factor*max(flick_xpiv(:,1));
    xLmax = parameters.piv_data_filename_interior.conversion_factor*min(flick_xpiv(:,end));
    yLmin = parameters.piv_data_filename_interior.conversion_factor*max(flick_ypiv(1,:)); 
    yLmax = parameters.piv_data_filename_interior.conversion_factor*min(flick_ypiv(end,:)); 
else
    xLmin = parameters.domainlimits(1); 
    xLmax = parameters.domainlimits(2); 
    yLmin = parameters.domainlimits(3); 
    yLmax = parameters.domainlimits(4); 
end 

%setting x_length and y_length 

xlength = xLmax - xLmin; 
ylength = yLmax - yLmin; 


coarsedx = xlength/parameters.Nxcoarse;
coarsedy = coarsedx;

parameters.Nycoarse = floor(ylength/coarsedy);
parameters.Nycoarse_shift = coarsedy*(ylength/coarsedy - floor(ylength/coarsedy))/2;

parameters.dx = xlength/parameters.Nx;
parameters.dy = parameters.dx;

%minus 4 allows for extrap for weno2
parameters.xlength = (parameters.Nxcoarse-4)*coarsedx;
parameters.ylength = (parameters.Nycoarse-4)*coarsedy;

%flickdata shift 
parameters.xshift_piv_data(1) = -xLmin-2*coarsedx;
parameters.yshift_piv_data(1) = -yLmin-2*coarsedy-parameters.Nycoarse_shift;




