function save_data_hairs_c(paths, parameters, simulation, initial)
%function save_data_hairs_c.m 
%input: initial - if this is the first time data is being saved 
%output: 
%saves hairs_c data

%global run_id pcount pathbase_results
%global ptindex_hairs hairs_c hairs_center


T1 = evalc(strcat('hairs_c_',num2str(simulation.pcount),' = simulation.hairs_c'));
T2 = evalc(strcat('hairs_a_',num2str(simulation.pcount),' = simulation.hairs_a'));
%ptindex_hairs = simulation.ptindex_hairs;
%hairs_center = simulation.hairs_center;

%the default extension is .mat
filename1 = [paths.pathbase_results, 'hairs_c_', parameters.run_id];
filename2 = [paths.pathbase_results, 'hairs_a_', parameters.run_id];

if (initial == 1) %if this is the first time saving then create the file
  save(filename1, 'hairs_c_*');
  save(filename2, 'hairs_a_*');
elseif (initial == 0) %if not the first time then append to the file
   save(filename1, 'hairs_c_*', '-append');
   save(filename2, 'hairs_a_*', '-append');
end