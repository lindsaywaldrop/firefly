function save_data_c(paths, parameters, simulation, initial)
%function save_data_c.m 
%input: initial - if this is the first time data is being saved 
%output: 
%saves c data

%global run_id pcount pathbase_results
%global c


T = evalc(['c_' num2str(simulation.pcount) ' = simulation.c']);

%the default extension is .mat
filename = [paths.pathbase_results,'c_', parameters.run_id];

if (initial == 1) %if this is the first time saving then create the file
  save(filename, 'c_*');
elseif (initial == 0) %if not the first time then append to the file
   save(filename, 'c_*', '-append');
end


