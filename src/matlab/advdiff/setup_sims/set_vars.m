%script set_vars.m
%initializes the variables for the program

%global pathbase_data hairNum fluid pathbase_piv piv_data_filename
%global xlength ylength Nx Ny dx dy x y u v
%global dthairfactor dt dt_flick t_final_flick t_steps_flick t_steps
%global dt dt_flick t_final_flick t_steps_flick t_steps
%global list_print_times t pcount print_time run_id
%global D initc domainlimits
%global c c_diff_matrix_l c_diff_matrix_u c_diff_RHS c_diff_matrix usegmres
%global explicit_vel piv_data_filename piv_data_filename_interior 
%global xshift_piv_data yshift_piv_data
%global handle_hairs hairs_data_filename hairs_data_filename_interior 
%only used in read_in_velocity_data_p2.m and advect_c.m
%global uplusx_piv uminusx_piv uplusy_piv uminusy_piv
%global vplusx_piv vminusx_piv vplusy_piv vminusy_piv
%global uplus2x_piv uminus2x_piv vplus2ypiv vminus2ypiv
%only used in read_in_velocity_data_p1.m and read_in_velocity_data_p2.m
%global Nxcoarse
%only used in advect_c if using weno
%global weno_eps
%only used in dealing with the hairs in setup_hairs.m and
%conecentration_absorbed_by_hairs
%global ptindex_hairs hairs_c hairs_center allhairs_center shift_hairs far_right_hair
%c bc on the right wall if using dirichlet bcs - always used for advection
%and for diffusion if specified 
%global cplusx_dbc diffusionrhsbc_flick

%setting the path so matlab can find all the functions
% addpath('./save_functions')
 
%read in the parameters
parameters.run_id = filenumber;
simulation.run_id = parameters.run_id;

run odorcapture_params.m

%cd(strcat(pathbase_data,'hairinfo-files/',num2str(hairNum),'hair_files/'))
%set x_length and y_length here if based on experimental data otherwise in
%params file manually
if strcmp(parameters.explicit_vel,'piv_data')
   if (parameters.handle_hairs)
       disp('setup_hairs')
       [parameters, simulation] = setup_hairs_for_velocity(paths, parameters);
   end
   disp('read_in_velocities')
   [parameters, simulation] = read_in_velocity_data_p1(paths, parameters, simulation);
end

%initializes Nx and Ny 
parameters.Nx = round(parameters.xlength/parameters.dx);
parameters.Ny = round(parameters.ylength/parameters.dy);

%should already be an int from read_in_velocity_data_p1.m
if (abs(parameters.Nx*parameters.dx-parameters.xlength)>1e-15) || (abs(parameters.Ny*parameters.dy-parameters.ylength)>1e-15)
    error('domain setting in set_vars.m')
end

%initialize x and y 
simulation.x(1:parameters.Nx+1,1)=(0:parameters.Nx)*parameters.dx; 
simulation.y(1:parameters.Ny+1,1)=(0:parameters.Ny)*parameters.dy; 
%for periodic/noflux bc 
%x(1:Nx,1)=(0:Nx-1)*dx; 
%y(1:Ny+1,1)=(0:Ny)*dy; 

%initializing time stepping variables
parameters.dt = min(0.9*parameters.dx/parameters.dtfactor, parameters.dthairfactor^2/4/parameters.D); 
if exist('parameters.dtmultiplier')
    parameters.dt = parameters.dtmultiplier*parameters.dt; 
end
%dt_rest = min(0.9*dx,dthairfactor^2/4/D); 

%flick time
if exist('parameters.t_final_factor_flick', 'var')
   parameters.t_final_flick = parameters.dt*parameters.t_final_factor_flick;
end
parameters.t_steps_flick = ceil(parameters.t_final_flick/parameters.dt);
parameters.dt_flick = parameters.t_final_flick/parameters.t_steps_flick;

simulation.t = 0; 
simulation.pcount = 1; 
simulation.list_print_times = 0;  
simulation.t_steps = 0; 

%initialize u and v 
simulation.u = zeros(parameters.Nx+1,parameters.Ny+1);
simulation.v = zeros(parameters.Nx+1,parameters.Ny+1); 
%for periodic/noflux bc
%u = zeros(Nx,Ny+1);
%v = zeros(Nx,Ny+1); 

%printing informaton about grids
%fprintf('\t dx = %4.16f\n\t dy = %4.16f\n\t dt = %4.25f\n\t dt_rest = %4.25f\n ', dx, dy, dt, dt_rest);
%fprintf('\t Grid (Nx by Ny) : %d by %d\n', Nx, Ny);
%fprintf('\t final times: flick: %4.16f\t, t_final_flick);
%fprintf('\t number of tsteps: flick: %d\t,t_steps_flick);
%fprintf('\t xlength by ylength: %4.16f by %4.16f\n',xlength,ylength); 

fprintf('\t dx = %4.16f\n\t dy = %4.16f\n\t dt = %4.25f\n ', parameters.dx, parameters.dy, parameters.dt_flick);
fprintf('\t Grid (Nx by Ny) : %d by %d\n', parameters.Nx, parameters.Ny);
fprintf('\t final times: flick: %4.16f\n', parameters.t_final_flick);
fprintf('\t number of tsteps: flick: %d\n',parameters.t_steps_flick);
fprintf('\t xlength by ylength: %4.16f by %4.16f\n',parameters.xlength,parameters.ylength); 

if (parameters.handle_hairs)
    [parameters, simulation] = setup_hairs(paths, parameters, simulation); 
end


%initialize the concentration
[parameters, simulation] = initialize_c(paths, parameters, simulation);

