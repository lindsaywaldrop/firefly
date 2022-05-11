function save_data_vel(initial, flickorreturn, paths, parameters, velocities)
%function save_data_vel.m 
%input: 
%output: 
%saves velocity data

%global run_id pathbase_results
%global u v 

%the default extension is .mat

filename = strcat(paths.pathbase_results, 'velocity_', parameters.run_id);

T = evalc(['u_' flickorreturn ' = velocities.u']);
T2 = evalc(['v_' flickorreturn ' = velocities.v']);


if (initial == 1)
    save(filename, 'u_*','v_*');
elseif (initial == 0)
    save(filename, 'u_*', 'v_*', '-append');
end
      


