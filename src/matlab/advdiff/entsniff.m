function entsniff(topdir, Species, filenumbers, fluid, clpool)
% entsniff.m
%
% Script for running code on Bridges
%
%	This script imports SAMRAI data into Matlab, interpolates to a
%	uniform grid. Only the final time step (set by user) is imported.
%
%   Note: This script relies on Eric Tytell's Matlab code for importing samrai data! 
%	This package needs to be installed for this to work properly!'
% 
% 	GridSize = IBAMR fluid grid size
%	final_time = final time step of simulation, or time step of interest.
% 	n = total number of simulations
%	pathbase = where runs are located

paths.topdir = topdir;
%Change directory
cd(paths.topdir)

% Add paths to relevant matlab analysis scripts
addpath(genpath(strcat(topdir,'/src/matlab/advdiff')))

% global pathbase_piv pathbase_data pathbase_results GridSize final_time 
% global files files0 hairNum fluid topdir

parameters.GridSize = 512;
parameters.final_time = 5000;
parameters.fluid = fluid;  % options: 'water' or 'air'
%parameters.D = 6.02e-7;       %current: 6.02e-5 in air; prev: 6.02e-2caproic acid in air  - in m^2/s -> correspond to half_exp IC 
parameters.D = 1e-11;    

j = length(filenumbers);
for ii = 1:j
   files{ii} = sprintf('%d', filenumbers(ii));
   files0{ii} = sprintf('%04d', filenumbers(ii));
end

% Setting paths to necessary files
paths.pathbase_ibamr = strcat(topdir, '/results/ibamr/', Species, '/');
paths.pathbase_data = strcat(topdir, '/data/');
paths.pathbase_results = strcat(topdir, '/results/odorcapture/', Species, '/');
parameters.Species = Species;
parameters.L = 400e-6;

if clpool == 1
    
    for i=1:length(files)
        % i=3
        % tic
        disp(['Simulation number: ', files{i}])
        disp('   ')
        
 		% Setting up hair info files
        if isfile(strcat(paths.pathbase_data, 'hairinfo-files/', parameters.Species,...
			 	'/hairinfo', files{i}, '.mat')) == 0
			 delete(strcat(paths.pathbase_data, 'hairinfo-files/', parameters.Species,...
			 		'/hairinfo', files{i}, '.mat'))
        end
        
        disp(['Setting up hair info files for ', files{i}])
		convert_hairdata(paths.pathbase_data, parameters.Species, str2double(files{i}));
    
        if isfile(strcat(paths.pathbase_ibamr, 'viz_IB2d', files{i}, '.mat')) == 0
            delete strcat(paths.pathbase_ibamr, 'viz_IB2d', files{i}, '.mat')
        end
        
		% Interpolates velocity fields and saves.
        disp(['Interpolating velocity fields for ', files{i}])
        entsniffinterp(i, files, paths.pathbase_ibamr, parameters.GridSize, sprintf('%05d',parameters.final_time));
        disp('    ')

        % Run Shilpa's code
        disp(['starting simulation for ', files0{i}])
        crabs(paths, parameters, files0{i})
        
    end
    
elseif clpool > 1
    
    mycluster = parpool(clpool);
    
    parfor i = 1:length(files)
        % tic
        disp(['Simulation number: ', files{i}])
        disp('   ')
        
 		% Setting up hair info files
        if isfile(strcat(paths.pathbase_data, 'hairinfo-files/', parameters.Species,...
			 	'/hairinfo', files{i}, '.mat')) == 0
			 delete(strcat(paths.pathbase_data, 'hairinfo-files/', parameters.Species,...
			 		'/hairinfo', files{i}, '.mat'))
        end
        
        disp(['Setting up hair info files for ', files{i}])
		convert_hairdata(paths.pathbase_data, parameters.Species, str2double(files{i}));
    
        if isfile(strcat(paths.pathbase_ibamr, 'viz_IB2d', files{i}, '.mat')) == 0
            delete strcat(paths.pathbase_ibamr, 'viz_IB2d', files{i}, '.mat')
        end
        
		% Interpolates velocity fields and saves.
        disp(['Interpolating velocity fields for ', files{i}])
        entsniffinterp(i, files, paths.pathbase_ibamr, parameters.GridSize, sprintf('%05d',parameters.final_time));
        disp('    ')

        % Run Shilpa's code
        disp(['starting simulation for ', files0{i}])
        crabs(paths, parameters, files0{i})
        
    end
else 
    disp('Error: Not a valid value for clpool!')
end

end
