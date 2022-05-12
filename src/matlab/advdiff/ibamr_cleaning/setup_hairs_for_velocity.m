function [parameters] = setup_hairs_for_velocity(paths, parameters) 

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

run_id = str2double(parameters.run_id);
all_pts_file = strcat(paths.pathbase_data,'vertex-files/',parameters.Species,'/',...
    parameters.Species,'_',num2str(run_id),'.vertex');

all_pts = dlmread(all_pts_file, ' ');
all_pts = all_pts(2:end,:);
min_x = min(all_pts(:,1));
max_x = max(all_pts(:,1));

max_x_hair = zeros(1,parameters.hairNum);

for peep=1:parameters.hairNum
    pts_x = eval(['p.h',num2str(peep),'_x']);
    pts_y = eval(['p.h',num2str(peep),'_y']);
    max_x_hair(peep) = max(max(pts_x));
    blip = boundary(pts_x,pts_y);
    h1_x = pts_x(blip);
    h1_y = pts_y(blip);
    flick_x_hairs{peep} = parameters.hairs_data_filename_interior.conversion_factor.*h1_x;
    flick_y_hairs{peep} = parameters.hairs_data_filename_interior.conversion_factor.*h1_y;
end

parameters.flick_x_hairs = flick_x_hairs;
parameters.flick_y_hairs = flick_y_hairs;

parameters.far_right_hair = max(max_x_hair);

parameters.dthairfactor = 1; 
%loads the location of the hairs
%flick_hairs = eval(['flickdata.' parameters.hairs_data_filename_interior.hairs]);

% Set domain limits based on hairs
%xLmin = min_x - 0.05*parameters.L;
%xLmax = max_x + 0.15*parameters.L;
%yLmin = 10e-6;
%yLmax = 0.5*parameters.L;
%yLmin = min([flick_hairs.y]) - 0.05;
%yLmax = max([flick_hairs.y]) + 0.05;
%if (xLmax > 0.5*parameters.L) 
%    xLmax = 0.5*parameters.L;
%end
%if (xLmin < -0.5*parameters.L)
%    xLmin = -0.5*parameters.L;
%end
%parameters.domainlimits = [xLmin, xLmax, yLmin, yLmax];

% if (parameters.hairs_data_filename_interior.givenradius)
%     
%      for i=1:length(flick_hairs)
%         %flick hairs centerlocation
%         simulation.flick_x_hairs = parameters.hairs_data_filename_interior.conversion_factor*flick_hairs(i).x;
%         simulation.flick_y_hairs = parameters.hairs_data_filename_interior.conversion_factor*flick_hairs(i).y;
%         simulation.flick_hairs_center(i,:) = [simulation.flick_x_hairs, simulation.flick_y_hairs];
%         
%      end
%      
%  
%       parameters.allhairs_center(1,:) = mean(simulation.flick_hairs_center);
%       parameters.dthairfactor = eval(['flickdata.' parameters.hairs_data_filename_interior.radius]); 
%         
%  end