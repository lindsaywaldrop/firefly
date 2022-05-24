function [parameters, simulation] = initialize_c(paths, parameters, simulation)
%function initialize_c.m 
%input: initc - string stating which initialization is desired
%output: 
%sets c initially

%global c x y xlength ylength
%global far_right_hair
%global cplusx_dbc diffusionrhsbc_flick

[xx,yy] = ndgrid(simulation.x, simulation.y);


if strcmp(parameters.initc,'constant')
    %constant
    simulation.c = xx*0+0.5;
    
elseif strcmp(parameters.initc,'exp_right_small') %THIS ONE 
    %NOTE THAT THIS IS NOW THE SAME FOR SMALL AND LARGE CASES 
    %width = 1.2;
    %c_Linf = 7; 
    %exp_center = 1.45;
    %c_max = 0.25; 
    
    width = 0.050*parameters.L; 
    runid = str2double(simulation.run_id);
    hair_vertices = dlmread(strcat(paths.pathbase_data,'vertex-files/',parameters.Species,...
        '/',parameters.Species,'_', num2str(runid),'.vertex'),' ');
    far_right = max(max(hair_vertices(2:end,1)))
    %%dist_frh = 0.0125;
    %dist_frh = 0.005;
    
    %if parameters.far_right_hair <= 0.005990
    %	exp_center = 0.005990+width/2; %far_right_hair + dist_frh + width/2; 
    %else
	exp_center = far_right+width/2; %far_right_hair + dist_frh + width/2; 
	%end
    c_Linf = 7; 
    c_max_constant = 0.1; 
    parameters.c_max = c_max_constant/parameters.ylength; 
    
    simulation.c = ((xx >= exp_center-width/2)&(xx <= exp_center+width/2)).*parameters.c_max.*exp((-c_Linf*((2*(xx-exp_center)/width)).^2));
    parameters.c_max
    max(max(simulation.c))
    
    parameters.cplusx_dbc = 0; 
    parameters.diffusionrhsbc_flick = 'noflux'; 
    
elseif strcmp(parameters.initc,'half_exp') %THIS ONE 
    width = 0.05*parameters.L; 
    runid = str2double(simulation.run_id);
    hair_vertices = dlmread(strcat(paths.pathbase_data,'vertex-files/',parameters.Species,...
        '/',parameters.Species,'_', num2str(runid),'.vertex'),' ');
    far_right = max(max(hair_vertices(2:end,1)))
    %%dist_frh = 0.0125;
    %dist_frh = 0.005;
    %if parameters.far_right_hair <= 0.005990
    %	exp_center = 0.005990+width/2; %far_right_hair + dist_frh + width/2; 
    %else
	exp_center = far_right+width*2 %far_right_hair + dist_frh + width/2; 
	%end
    c_Linf = 7; 
    c_max_constant = 0.1; 
    parameters.c_max = c_max_constant/parameters.ylength; 
    
    simulation.c = (((xx >= exp_center-width/2)&(xx <= exp_center)).*parameters.c_max.*exp((-c_Linf*((2*(xx-exp_center)/width)).^2))) + (xx>exp_center).*parameters.c_max;    
    parameters.c_max
    max(max(simulation.c))
    
    parameters.cplusx_dbc = parameters.c_max; 
    parameters.diffusionrhsbc_flick = 'dirichlet'; 
    
else
    error('initc is not a valid choice!');
end
 
%c = cos(xx*pi) + 1;     
%C = exp(-(2*pi*xx-pi).^2) + exp(-(2*pi*yy-pi).^2); 
%C = exp(-(xx-pi).^2) + exp(-(yy-pi).^2); 
%C = exp(-(xx-pi).^2).*(1+cos(yy)); 
%C = cos(yy) + 1; 
%C = exp(-(xx-pi).^2);
