function save_printdata(paths, parameters, simulation)
%function save_printdata.m 
%input:  
%output: 
%saves printing information data

pcount = simulation.pcount;
list_print_times = simulation.list_print_times;
t_steps = simulation.t_steps;
print_time = parameters.print_time;

%the default extension is .mat
filename = [paths.pathbase_results, 'initdata_', parameters.run_id];
save (filename, 'pcount', 'list_print_times', 't_steps', 'print_time', '-append');
