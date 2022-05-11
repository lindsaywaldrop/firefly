%PARAMETERS


%Nx - set grid only for x direction and it will define dx=dy and Ny 
   Nx = 512;

%for piv data as velocity field
%only used in read_in_velocity_data_p1.m 
%Nxcoarse must be a factor of Nx - see read_in_velocity_data_p1.m
Nxcoarse = Nx; %right now set the same as Nx 

%time stepping factor dt = 0.9*dx/dtfactor - need to set this using the max
%velocity
   dtfactor = 0.15;   %0.15 m/s -> max velocity in all the exp data
   %dthairfactor = 0.0045; %hair length scale in m - now set when loading hairinfo

   %final time (s):                 
   t_final_flick = 20; %0.1 s -> 200 s 
   %t_final_factor_flick = 20000;

%print every print_time timesteps 
print_time = 10;  
 
%if using weno for advection then need to set such that 
%weno_eps = 1e-6*O(u^2)
   weno_eps = 1e-6;

%initializing the bulk surfactants
%initc = 'exp_right_small';
initc = 'half_exp';

%parameters
explicit_vel = 'piv_data';  
                            
%domain size 
%if want automatically set from piv_data then set doaminlimits = 'auto'
  %otherwise domainlimits = [xLmin, xLmax, yLmin, yLmax] in the correctly
%converted units 
  domainlimits = [-0.55, 0.7, -0.625, 0.625]; %in m 

%only needed if explicit_vel = 'piv_data'                            
  piv_data_filename = 'viz_IB2d1233';
%piv_data_returnfilename = 'data_being_used/simdata/marine_water/set1-return/marine-water_returndata_shilpa';
%this information needs to hold for both data files
%piv_data_filename_interior.filename = 'flickdata'; 
piv_data_filename_interior.x = 'x';
piv_data_filename_interior.y = 'y';
piv_data_filename_interior.u = 'Vinterp.U_x';
piv_data_filename_interior.v = 'Vinterp.U_y';
piv_data_filename_interior.conversion_factor = 1; 
piv_data_filename_interior.forcedivfree = 0; 

%diffusion solver 
usegmres = 0; %0 if LU decomposition and 1 if gmres iterative solver
                            
%diffusion coefficient (m^2/s)
%D = 2000*7.84e-10;     %caproic acid in water - in m^2/s m-> corresponds to exp_right_small IC 
D = 2000*6.02e-6;       %caproic acid in air  - in m^2/s -> correspond to half_exp IC 
%D = 5e-4;        %testing

%hairs  
handle_hairs = 1; 
hairs_data_filename = 'hairinfo1233';
%hairs_data_returnfilename = 'data_being_used/simdata/marine_water/set1-return/marine-water_returndata_shilpa';
%this information needs to hold for both data files
hairs_data_filename_interior.filename = 'p';
hairs_data_filename_interior.numofhairs = 3; 
hairs_data_filename_interior.hairs = 'hairs'; 
hairs_data_filename_interior.givenradius = 1; 
hairs_data_filename_interior.radius = 'radius_m'; 
hairs_data_filename_interior.conversion_factor = 1; 
