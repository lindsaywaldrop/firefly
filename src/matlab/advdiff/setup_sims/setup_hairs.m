function [parameters, simulation] = setup_hairs(paths, parameters, simulation) 

%global hairs_data_filename hairs_data_filename_interior
%global pathbase_data hairNum run_id
%global xshift_piv_data yshift_piv_data
%global x y domainlimits
%global ptindex_hairs hairs_c hairs_center allhairs_center shift_hairs
%global far_right_hair   

%NOTES: SK 2018_01_22
%COMMENTED ALL RETURN HAIR STUFF OUT 
%AND THEREFORE ALSO NO NEED TO SHIFT SINCE ONLY FLICK 
%ASSUMED GIVEN RADIUS OF HAIRS 

[xx,yy] = ndgrid(simulation.x, simulation.y);

%hairs filename
flickdata = load(strcat(paths.pathbase_data, 'hairinfo-files/', num2str(parameters.hairNum),...
    'hair_files/', parameters.hairs_data_filename, '.mat'));  
% [flickdata.p] = convert_hairdata(pathbase_data,hairNum,str2double(run_id));


%loads the location of the hairs
flick_hairs = eval(['flickdata.' parameters.hairs_data_filename_interior.filename '.' parameters.hairs_data_filename_interior.hairs]); 

% Set domain limits based on hairs
xLmin = min([flick_hairs.x]) - 0.05;
xLmax = max([flick_hairs.x]) + 0.15;
yLmin = min([flick_hairs.y]) - 0.05;
yLmax = max([flick_hairs.y]) + 0.05;
parameters.domainlimits = [xLmin, xLmax, yLmin, yLmax];

%initializes the concentration absorbed by each hair
simulation.hairs_c = struct([]);

%if a radius is given for each hairs and the location of the hairs just
%gives the center of each hair 
if (parameters.hairs_data_filename_interior.givenradius)
   
    %%NO LONGER TRUE  
    %%we assume here that the hairs in the flick and return have the same size and shape
    %%but the code is written such that the location of the center of all
    %%hairs could be different 
    
    %if (length(flick_hairs)~=length(return_hairs))
    %    error('number of hairs is different in flick and return') 
    %end

    %we use the information from the flick hairs 
    %radius of the hairs from params file -> now in meters 
     
    parameters.hairs_radius = parameters.hairs_data_filename_interior.conversion_factor*...
    	eval(['flickdata.'  parameters.hairs_data_filename_interior.filename '.' parameters.hairs_data_filename_interior.radius]);
    %hairs_radius = 4*hairs_radius
    
    %find indices of points inside the hairs 
    for i=1:length(flick_hairs)
        distance_to_hair_squared = 0;
        simulation.x_hairs = 0;
        simulation.y_hairs = 0; 
        simulation.x_hairs = parameters.hairs_data_filename_interior.conversion_factor*flick_hairs(i).x + parameters.xshift_piv_data(1);
        simulation.y_hairs = parameters.hairs_data_filename_interior.conversion_factor*flick_hairs(i).y + parameters.yshift_piv_data(1); 
		distance_to_hair_squared = (xx-simulation.x_hairs).^2 + (yy-simulation.y_hairs).^2;
        %linear indices in matrix c 
        [simulation.ptindex_hairs{i}]= find(distance_to_hair_squared <= parameters.hairs_radius^2); 
        %finishes initialization of the concentration absorbed by each hair
        simulation.hairs_c{i} = zeros(length(simulation.ptindex_hairs{i}),1); 
        simulation.hairs_center(i,:) = [simulation.x_hairs, simulation.y_hairs];
    end
    %finding the farthest right point on any hair based on just the flick hairs
    parameters.far_right_hair = max(simulation.hairs_center(:,1))+parameters.hairs_radius; 
    
%    
    
    
    
end
    


% figure(1)
% hold on
% plot(xx,yy,'bx')
% for i=1:length(flick_hairs)
%     plot(xx(ptindex_hairs{i}),yy(ptindex_hairs{i}),'ro')
%     plot(hairs_center(i,1),hairs_center(i,2),'gx')
%     if (hairs_data_filename_interior.givenradius)
%         viscircles(hairs_center(i,:),hairs_radius,'EdgeColor','k','LineWidth',1)
%     elseif (~hairs_data_filename_interior.givenradius)
%         plot(hairs_data_filename_interior.conversion_factor*hairs(i).x + xshift_piv_data,hairs_data_filename_interior.conversion_factor*hairs(i).y + yshift_piv_data,'k')
%         plot(hairs_data_filename_interior.conversion_factor*hairs(i).x + xshift_piv_data,hairs_data_filename_interior.conversion_factor*hairs(i).y + yshift_piv_data,'k*')
%     end
% end
% pause

    
