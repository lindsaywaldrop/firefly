%PARAMETERS


%Nx - set grid only for x direction and it will define dx=dy and Ny 
parameters.Nx = parameters.GridSize;
parameters.fluid = 'air';

%for piv data as velocity field
%only used in read_in_velocity_data_p1.m 
%Nxcoarse must be a factor of Nx - see read_in_velocity_data_p1.m
parameters.Nxcoarse = parameters.Nx; %right now set the same as Nx 

%time stepping factor dt = 0.9*dx/dtfactor - need to set this using the max
%velocity
parameters.dtfactor = 1;   %0.15 m/s -> max velocity in all the exp data
%dthairfactor = 0.0045; %hair length scale in m - now set when loading hairinfo

%if using weno for advection then need to set such that 
%weno_eps = 1e-6*O(u^2)
parameters.weno_eps = 1e-6;

%parameters
parameters.explicit_vel = 'ibamr_data';  
                            
%domain size 
%if want automatically set from piv_data then set doaminlimits = 'auto'
  %otherwise domainlimits = [xLmin, xLmax, yLmin, yLmax] in the correctly
%converted units 
% Note: domain limits are now defined in setup_hairs() and
% setup_hairs_for_velocity()
%domainlimits = [0, 0.3, -0.15, 0.22]; %in m original
%domainlimits = 'auto'; 
%only needed if explicit_vel = 'piv_data'                            
parameters.ibamr_data_filename = strcat('viz_IB2d', num2str(str2double(parameters.run_id)));
%piv_data_returnfilename = 'data_being_used/simdata/marine_water/set1-return/marine-water_returndata_shilpa';
%this information needs to hold for both data files
%piv_data_filename_interior.filename = 'flickdata'; 
parameters.ibamr_data_filename_interior.x = 'x';
parameters.ibamr_data_filename_interior.y = 'y';
parameters.ibamr_data_filename_interior.u = 'Vinterp.U_x';
parameters.ibamr_data_filename_interior.v = 'Vinterp.U_y';
parameters.ibamr_data_filename_interior.conversion_factor = 1; 
parameters.ibamr_data_filename_interior.forcedivfree = 0; 

%diffusion solver 
parameters.usegmres = 0; %0 if LU decomposition and 1 if gmres iterative solver
                     

if strcmp(parameters.fluid,'air')
	%initializing the bulk surfactants
	parameters.initc = 'half_exp';
	%diffusion coefficient (m^2/s)	
	parameters.D = 2000*6.02e-6;       %caproic acid in air  - in m^2/s -> correspond to half_exp IC 
	simulation.D = parameters.D;
	%print every print_time timesteps 
	parameters.print_time = 500;  
	simulation.print_time = parameters.print_time;
	%final time (s):                 
	parameters.t_final_flick = 0.5; %0.1 s -> 200 s 
	simulation.t_final_flick = parameters.t_final_flick;
	%t_final_factor_flick = 20000;

elseif strcmp(parameters.fluid,'water')
	%initializing the bulk surfactants
	parameters.initc = 'exp_right_small';
	%diffusion coefficient (m^2/s)
	parameters.D = 2000*7.84e-10;     %caproic acid in water - in m^2/s -> corresponds to exp_right_small IC 
	simulation.D = parameters.D;
	%print every print_time timesteps 
	parameters.print_time = 100;  
	simulation.print_time = parameters.print_time;
	%final time (s):                 
	parameters.t_final_flick = 15; %0.1 s -> 200 s 
	simulation.t_final_flick = parameters.t_final_flick;
   %t_final_factor_flick = 20000;
else 
	disp('unknown fluid type')
	return
end

disp(['Running in ',parameters.fluid,' with D=', num2str(parameters.D), ' for ', num2str(parameters.t_final_flick),' s.'])
disp(' ')

%hairs  
parameters.handle_hairs = 1; 
parameters.hairs_data_filename = strcat('hairinfo',num2str(str2double(parameters.run_id)));
%this information needs to hold for both data files
parameters.hairs_data_filename_interior.filename = 'p';
parameters.hairs_data_filename_interior.numofhairs = parameters.hairNum; 
parameters.hairs_data_filename_interior.hairs = 'hairs'; 
%parameters.hairs_data_filename_interior.givenradius = 0.001; 
%parameters.hairs_data_filename_interior.radius = 'radius_m'; 
parameters.hairs_data_filename_interior.conversion_factor = 1; 
