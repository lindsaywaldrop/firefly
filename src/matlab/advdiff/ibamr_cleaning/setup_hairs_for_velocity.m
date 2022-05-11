function [parameters, simulation] = setup_hairs_for_velocity(paths, parameters) 

%global hairs_data_filename hairs_data_filename_interior 
%global allhairs_center pathbase_piv pathbase_data
%global dthairfactor run_id hairNum domainlimits

%NOTES: SK 2017_11_16
%COMMENTED ALL RETURN HAIR STUFF OUT 
%AND THEREFORE ALSO NO NEED TO SHIFT SINCE ONLY FLICK 
%ASSUMED GIVEN RADIUS OF HAIRS 

%hairs filename
flickdata_temp = load(strcat(paths.pathbase_data, '/hairinfo-files/', parameters.Species, '/', ...
	parameters.hairs_data_filename, '.mat'));  
flickdata = eval(['flickdata_temp.' parameters.hairs_data_filename_interior.filename]);

%flickdata = modify_hair_data(flickdata,parameters.hairs_data_filename_interior.numofhairs); 
eval([parameters.hairs_data_filename_interior.filename, '=flickdata']); 
save(strcat(paths.pathbase_data, '/hairinfo-files/', parameters.Species, '/', ...
	parameters.hairs_data_filename, '.mat'), parameters.hairs_data_filename_interior.filename); 

%loads the location of the hairs
flick_hairs = eval(['flickdata.' parameters.hairs_data_filename_interior.hairs])

% Set domain limits based on hairs
%xLmin = min([flick_hairs.x]) - 0.05;
%xLmax = max([flick_hairs.x]) + 0.15;
%yLmin = min([flick_hairs.y]) - 0.05;
%yLmax = max([flick_hairs.y]) + 0.05;
%parameters.domainlimits = [xLmin, xLmax, yLmin, yLmax];

if (parameters.hairs_data_filename_interior.givenradius)
    
     for i=1:length(flick_hairs)
        %flick hairs centerlocation
        simulation.flick_x_hairs = parameters.hairs_data_filename_interior.conversion_factor*flick_hairs(i).x;
        simulation.flick_y_hairs = parameters.hairs_data_filename_interior.conversion_factor*flick_hairs(i).y;
        simulation.flick_hairs_center(i,:) = [simulation.flick_x_hairs, simulation.flick_y_hairs];
        
     end
     
 
      parameters.allhairs_center(1,:) = mean(simulation.flick_hairs_center);
      parameters.dthairfactor = eval(['flickdata.' parameters.hairs_data_filename_interior.radius]); 
        
 end