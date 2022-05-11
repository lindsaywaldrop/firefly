function save_initdata(paths, parameters, simulation)
%function save_initdata.m 
%input: 
%output: 
%saves fixed parameters/variables initially

% Extracts values of usefulness.
dx = parameters.dx;
dy = parameters.dy;
Nx = parameters.Nx;
Ny = parameters.Ny; 
xlength = parameters.xlength; 
ylength = parameters.ylength; 
x = simulation.x; 
y = simulation.y; 
explicit_vel = parameters.explicit_vel;
D = parameters.D;
dt = parameters.dt; 
dt_flick = parameters.dt_flick;
t_final_flick = parameters.t_final_flick;
run_id = parameters.run_id;


%the default extension is .mat
filename = [paths.pathbase_results, 'initdata_', run_id, '.mat'];
save(filename, 'dx','dy', 'Nx', 'Ny', 'xlength', ...
     'ylength', 'x', 'y','run_id', 'explicit_vel', 'D',...
     'dt','dt_flick','t_final_flick');
