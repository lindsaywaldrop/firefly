function crabs(paths, parameters, filenumber)
%function crabs.m 
%input: filenumber - string
%output: 
%main file - calls all the timestepping steps for this program
%global topdir pathbase_data pathbase_piv pathbase_results hairNum
%global GridSize final_time fluid

%if isempty(paths.pathbase_data)
%    load(strcat(paths.topdir,'/src/matlab/temp_global_variable.mat'));
%end

%initialize
%simulation = struct([]);

% Setting parameters and simulation set up
set_vars 

%saving data initially
save_data(paths, parameters, simulation, 1); 

%FLICK
disp('Creating and saving velocity data')
disp('  ')
[velocities] = get_velocities(parameters.dt_flick/2,simulation.t,parameters.explicit_vel,'flick', paths, parameters, simulation);
save_data_vel(1, 'flick', paths, parameters, velocities);
disp('Done!')
disp(' ')

disp('Starting first flick...')

%time-stepping
%could speed this up by taking out the if else statements from the for loop

%advection - first step
[simulation] = advect_c(parameters.dt_flick/2, 'dirichlet', 'weno', ...
												parameters, simulation, velocities);
disp('.')
for timestep = 1:parameters.t_steps_flick
   
  %diffusion  
  %if first timestep then initialze the diffusion matrix   
  if (timestep == 1) 
      [simulation] = diffusion_c(parameters.dt_flick, 1, parameters.diffusionrhsbc_flick, parameters, simulation);
      disp('.')
  else
      [simulation] = diffusion_c(parameters.dt_flick, 0, parameters.diffusionrhsbc_flick, parameters, simulation); 
  end
  [simulation] = concentration_absorbed_by_hairs(simulation);
  %advection 
  %if not at the last timestep then step with dt but if at the last
  %timestep then step only dt/2    
  if (timestep ~= parameters.t_steps_flick)
     [simulation] = advect_c(parameters.dt_flick,'dirichlet','weno', ...
												parameters, simulation, velocities);
  elseif (timestep == parameters.t_steps_flick) 
     [simulation] = advect_c(parameters.dt_flick/2,'dirichlet','weno',  ...
												parameters, simulation, velocities);
  end
  
  simulation.t = simulation.t + parameters.dt_flick; 
  simulation.t_steps = simulation.t_steps + 1; 
  
  %saving data
  if (mod(simulation.t_steps, parameters.print_time)==0) 
    simulation.pcount = simulation.pcount + 1; 
    save_data(paths, parameters, simulation, 0);     
    simulation.list_print_times(simulation.pcount) = simulation.t; 
    fprintf('printing %g %g \n',simulation.t, simulation.pcount)
  end

end


%saving data finally
simulation.pcount = simulation.pcount+1;  
save_data(paths, parameters, simulation, 0);     
simulation.list_print_times(simulation.pcount) = simulation.t; 
fprintf('printing %g %g \n', simulation.t, simulation.pcount)
save_printdata(paths, parameters, simulation); 
 