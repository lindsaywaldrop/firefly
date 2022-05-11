function save_data_hairs_c(paths, parameters, simulation, initial)
%function save_data_hairs_c.m 
%input: initial - if this is the first time data is being saved 
%output: 
%saves hairs_c data

%global run_id pcount pathbase_results
%global ptindex_hairs hairs_c hairs_center


T = evalc(strcat('hairs_c_',num2str(simulation.pcount),' = simulation.hairs_c'));
ptindex_hairs = simulation.ptindex_hairs;
hairs_center = simulation.hairs_center;

%the default extension is .mat
filename = [paths.pathbase_results, 'hairs_c_', parameters.run_id];

if (initial == 1) %if this is the first time saving then create the file
  save(filename, 'hairs_c_*', 'ptindex_hairs','hairs_center');
elseif (initial == 0) %if not the first time then append to the file
   save(filename, 'hairs_c_*', '-append');
end