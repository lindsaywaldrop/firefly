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
elseif strcmp(parameters.initc,'constant_patch')
    %constant_patch - ISSUES WITH DISCONTINUITIES
    simulation.c = ((xx >= 1.25)&(xx <= 1.75))*0.5;
elseif strcmp(parameters.initc,'constant_patch_with_smoothing')
    %constant patch with smoothing  
    %exponential smoothing - NEEDS FIXING
    %x1 = 1.25; 
    %x2 = 1.75;
    %c_eps = 0.3;
    %c_Linf = 2;  %11
    %c_m = 0.5; 
    %c = ((xx >= x1)&(xx <= x2))*c_m + ...
    %    ((xx > x1-c_eps)&(xx < x1)).*c_m.*exp(-c_Linf*(x1-xx)/c_eps) + ...
    %    ((xx > x2)&(xx < x2 + c_eps)).*c_m.*exp(-c_Linf*(xx-x2)/c_eps); 
    x1 = 0.7; 
    x2 = 1.7;
    c_eps = 0.4;
    c_m = 0.5;
    xxl = (xx-x1)/c_eps;
    xxr = -(xx-x2)/c_eps;
    simulation.c = ((xx >= x1+c_eps)&(xx <= x2-c_eps))*c_m + ...
        ((xx > x1-c_eps)&(xx < x1+c_eps)).*c_m.*(1/2+(1/64)*(105*xxl - 175*xxl.^3 + 147*xxl.^5 - 45*xxl.^7)) + ...
        ((xx > x2-c_eps)&(xx < x2+c_eps)).*c_m.*(1/2+(1/64)*(105*xxr - 175*xxr.^3 + 147*xxr.^5 - 45*xxr.^7));
elseif strcmp(parameters.initc,'constant_patch_with_smoothing_tanh')
    %x1 = 0.7; 
    %x2 = 1.7;
    %c_eps = 0.4;
    x1 = 0.425; 
    x2 = 1.2;
    c_eps = 0.3875;
    %c_Linf = 2; %upwind advection convg test values
    %c_Linf = 7;  %diffusion convg test values 
    c_Linf = 4; 
    %c_m = 0.5;
    c_m = 0.25; 
    xxl = c_Linf*(xx-x1)/c_eps;
    xxr = -c_Linf*(xx-x2)/c_eps;
    %c = ((xx >= x1+c_eps)&(xx <= x2-c_eps))*c_m + ...
    %    ((xx > x1-c_eps)&(xx < x1+c_eps)).*c_m.*((1/2).*tanh(xxl)+(1/2))+...
    %    ((xx > x2-c_eps)&(xx < x2+c_eps)).*c_m.*((1/2).*tanh(xxr)+(1/2));
    simulation.c = ((xx >= x1-c_eps)&(xx <= x1+c_eps)).*c_m.*((1/2).*tanh(xxl)+(1/2))+...
        ((xx > x2-c_eps)&(xx <= x2+c_eps)).*c_m.*((1/2).*tanh(xxr)+(1/2));
elseif strcmp(parameters.initc,'cos_patch')  
    %cosine patch - ISSUES WITH DISCONTINUITIES
    simulation.c = ((xx >= 1.25)&(xx <= 1.75)).*((cos((xx-3/2)*4*pi)+1)*(1/4));
elseif strcmp(parameters.initc,'cos')
    simulation.c = (cos(xx*2*pi/parameters.xlength)+1)*(1/4);
elseif strcmp(parameters.initc,'exp')
    x1 = 0.2; 
    x2 = 1.4;
    %c_Linf = 4; %upwind advection convg test values
    c_Linf = 4; %weno advection convg test values
    %c_Linf = 10; %diffusion convg test values 
    %c = exp((-c_Linf*((2*(xx-x1)/(x2-x1)-1)).^2)); %upwind and diffusion convg test values
    simuation.c = 0.25*exp((-c_Linf*((2*(xx-xlength/2)/(x2-x1))).^2));
    %y1 = 0.2;
    %y2 = 1.2; 
    %c = 0.25*exp((-c_Linf*((2*(yy-ylength/2)/(y2-y1))).^2));
elseif strcmp(parameters.initc,'exp_right')
    x1 = 0.2; 
    x2 = 1.4;
    c_Linf = 4; 
    simulation.c = 1*exp((-c_Linf*((2*(xx-parameters.xlength/2-0.55)/(x2-x1))).^2));
elseif strcmp(parameters.initc,'exp_right_small') %THIS ONE 
    %NOTE THAT THIS IS NOW THE SAME FOR SMALL AND LARGE CASES 
    %width = 1.2;
    %c_Linf = 7; 
    %exp_center = 1.45;
    %c_max = 0.25; 
    
    width = 0.05; 
    %%dist_frh = 0.0125;
    %dist_frh = 0.005;
    
    if parameters.far_right_hair <= 0.005990
    	exp_center = 0.005990+width/2; %far_right_hair + dist_frh + width/2; 
    else
		exp_center = parameters.far_right_hair+width/2; %far_right_hair + dist_frh + width/2; 
	end
    c_Linf = 7; 
    c_max_constant = 0.1; 
    parameters.c_max = c_max_constant/parameters.ylength; 
    
    simulation.c = ((xx >= exp_center-width/2)&(xx <= exp_center+width/2)).*parameters.c_max.*exp((-c_Linf*((2*(xx-exp_center)/width)).^2));
    parameters.c_max
    max(max(simulation.c))
    
    parameters.cplusx_dbc = 0; 
    parameters.diffusionrhsbc_flick = 'noflux'; 
elseif strcmp(parameters.initc,'exp_right_small_v2') %THIS ONE 
    %NOTE THAT THIS IS NOW THE SAME FOR SMALL AND LARGE CASES 
    %width = 1.2;
    %c_Linf = 7; 
    %exp_center = 1.45;
    %c_max = 0.25; 
    
    width = 0.1; 
    xlength_air = 1.240234375000000;
    %%dist_frh = 0.0125;
    %dist_frh = 0.005;
    
    exp_center = parameters.far_right_hair+width/2; %far_right_hair + dist_frh + width/2; 
    c_Linf = 7; 
    c_max_constant = 0.1; 
    aa = sqrt(pi)*erf(sqrt(c_Linf))/4/sqrt(c_Linf);
    bb = xlength_air-exp_center+0.06*20;
    %c_max = (c_max_constant/ylength)*(bb/0.1/aa+1)/2/2;
    parameters.c_max = (c_max_constant/parameters.ylength)*(bb/0.1/aa+1)/2;
    
    simulation.c = ((xx >= exp_center-width/2)&(xx <= exp_center+width/2)).*parameters.c_max.*exp((-c_Linf*((2*(xx-exp_center)/width)).^2));
    parameters.c_max
    max(max(simulation.c))
    
    parameters.cplusx_dbc = 0; 
    parameters.diffusionrhsbc_flick = 'noflux'; 
    
elseif strcmp(parameters.initc,'exp_right_small_smdom') %THIS ONE 
    %NOTE THAT THIS IS NOW THE SAME FOR SMALL AND LARGE CASES 
    %width = 1.2;
    %c_Linf = 7; 
    %exp_center = 1.45;
    %c_max = 0.25; 
    
    width = 0.1; 
    %%dist_frh = 0.0125;
    %dist_frh = 0.005;
    
    exp_center = parameters.far_right_hair+width/2; %far_right_hair + dist_frh + width/2; 
    c_Linf = 7; 
    c_max_constant = 0.1; 
    parameters.c_max = c_max_constant/parameters.ylength; 
    
    simulation.c = ((xx >= exp_center-width/2)&(xx <= exp_center+width/2)).*parameters.c_max.*exp((-c_Linf*((2*(xx-exp_center)/width)).^2));
    parameters.c_max
    max(max(simulation.c))
    
    parameters.cplusx_dbc = 0; 
    parameters.diffusionrhsbc_flick = 'noflux'; 
    
elseif strcmp(parameters.initc,'exp_right_large') %THIS ONE - NOT USING THIS ONE  
    %width = 0.12;
    %c_Linf = 7; 
    %exp_center = 0.145;
    %c_max = 0.25; 
    
    width = 0.1; 
    dist_frh = 0.0125;
    exp_center = parameters.far_right_hair + dist_frh + width/2; 
    c_Linf = 7; 
    c_max_constant = 0.1; 
    parameters.c_max = c_max_constant/parameters.ylength; 
    
    simulation.c = ((xx >= exp_center-width/2)&(xx <= exp_center+width/2)).*parameters.c_max.*exp((-c_Linf*((2*(xx-exp_center)/width)).^2));  
    parameters.c_max
    max(max(simulation.c))
    
    parameters.cplusx_dbc = 0;
    parameters.diffusionrhsbc_flick = 'noflux'; 
    
elseif strcmp(parameters.initc,'exp_left')
    x1 = 0.2; 
    x2 = 1.4;
    c_Linf = 4; 
    c = 0.25*exp((-c_Linf*((2*(xx-parameters.xlength/2+0.15)/(x2-x1))).^2));
    
elseif strcmp(parameters.initc,'half_exp') %THIS ONE 
    width = 0.05; 
    %%dist_frh = 0.0125;
    %dist_frh = 0.005;
    if parameters.far_right_hair <= 0.005990
    	exp_center = 0.005990+width/2; %far_right_hair + dist_frh + width/2; 
    else
		exp_center = parameters.far_right_hair+width/2; %far_right_hair + dist_frh + width/2; 
	end
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
