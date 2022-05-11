function [newu,newv] = noslipboundary_givenradiushairs(bu,bv)

global hairs_data_filename hairs_data_filename_interior
global x y dx dy
global ptindex_hairs hairs_center
 
%we just use the information from the flick hairs 

regularized_fnc = ones(size(bu)); 
[xx,yy] = ndgrid(x,y);

num_hairs = length(ptindex_hairs); 

%hairs filename
load(['pivdata/' hairs_data_filename '.mat'])  
%radius of the hairs from params file -> given in cm and need in mm 
hairs_radius = hairs_data_filename_interior.conversion_factor*eval([hairs_data_filename_interior.filename '.' hairs_data_filename_interior.radius]);


%sets the velocity field to zero within each hair
for i=1:num_hairs
    regularized_fnc(ptindex_hairs{i}) = 0;  
end

%now regularizing to the outer velocity field
%heaviside_width*dx is the width of the regularization to 0 within the hair
%assumes dx = dy 
heaviside_width = 1; 

%find points external of hair but inside the regularized area
counter = 1; 
for i=1:num_hairs
    x_hairs = hairs_center(i,1);
    y_hairs = hairs_center(i,2); 
    distance_to_hair_squared = (xx-x_hairs).^2 + (yy-y_hairs).^2;
    [outerptindex_hairs{i}]= find((distance_to_hair_squared<=((hairs_radius+heaviside_width*dx)^2))&(distance_to_hair_squared>(hairs_radius^2))); 
    outerptindex_hairs_vector(counter:counter+length(outerptindex_hairs{i})-1,1) = outerptindex_hairs{i};
    outerptindex_whichhairs_vector(counter:counter+length(outerptindex_hairs{i})-1,1) = i;
    counter = counter+length(outerptindex_hairs{i});
end

%finding if two hairs have the same points in their regularized area 
[outerptindex_hairs_vector_sorted,outerptindex_hairs_vector_i] = sort(outerptindex_hairs_vector);

counter = 0;
for i=1:length(outerptindex_hairs_vector_sorted)-1
    if (outerptindex_hairs_vector_sorted(i) == outerptindex_hairs_vector_sorted(i+1))
        counter = counter + 1; 
        repeated_indices_sorted(counter) = i;
    end
end


if (counter > 0)
    repeated_indices = outerptindex_hairs_vector_i(repeated_indices_sorted);
    repeated_indices_plus = outerptindex_hairs_vector_i(repeated_indices_sorted+1);
    if (size(outerptindex_hairs_vector(repeated_indices)) ~= size(unique(outerptindex_hairs_vector(repeated_indices))))
        error('indices are repeated more than twice - noslipboundary_givenradiushairs.m')
    end

    %removing one of the repetitions
    for i=1:length(repeated_indices) 
    
        if (((xx(outerptindex_hairs_vector(repeated_indices(i)))-...
              hairs_center(outerptindex_whichhairs_vector(repeated_indices(i)),1))^2+...
            (yy(outerptindex_hairs_vector(repeated_indices(i)))-...
              hairs_center(outerptindex_whichhairs_vector(repeated_indices(i)),2))^2) <= ... 
           ((xx(outerptindex_hairs_vector(repeated_indices_plus(i)))-...
              hairs_center(outerptindex_whichhairs_vector(repeated_indices_plus(i)),1))^2+...
            (yy(outerptindex_hairs_vector(repeated_indices_plus(i)))-...
              hairs_center(outerptindex_whichhairs_vector(repeated_indices_plus(i)),2))^2))
    
            rmindex = find(outerptindex_hairs{outerptindex_whichhairs_vector(repeated_indices_plus(i))} ...
                           == outerptindex_hairs_vector(repeated_indices_plus(i)));
            outerptindex_hairs{outerptindex_whichhairs_vector(repeated_indices_plus(i))}(rmindex) = [];  
       
        else
            rmindex = find(outerptindex_hairs{outerptindex_whichhairs_vector(repeated_indices(i))} ...
                           == outerptindex_hairs_vector(repeated_indices(i)));                       
            outerptindex_hairs{outerptindex_whichhairs_vector(repeated_indices(i))}(rmindex) = [];  
        end
    end
end

for i=1:num_hairs
    x_hairs = hairs_center(i,1);
    y_hairs = hairs_center(i,2);
    for j=1:length(outerptindex_hairs{i})
        %finding the point on the outer circle that is normal from the center of
        %the hair
        vec_to_center_x = xx(outerptindex_hairs{i}(j))-x_hairs;
        vec_to_center_y = yy(outerptindex_hairs{i}(j))-y_hairs;
        dist_to_center = sqrt(vec_to_center_x^2+vec_to_center_y^2); 
        xp = ((hairs_radius+heaviside_width*dx)*vec_to_center_x/dist_to_center)+x_hairs;
        yp = ((hairs_radius+heaviside_width*dx)*vec_to_center_y/dist_to_center)+y_hairs; 
        
        vec_to_xp = xp-xx(outerptindex_hairs{i}(j));
        vec_to_yp = yp-yy(outerptindex_hairs{i}(j));
        dist_to_xpyp = sqrt(vec_to_xp^2+vec_to_yp^2);
        dist_to_hair = heaviside_width*dx - dist_to_xpyp; 
        if (dist_to_hair < 0)
            error('dist_to_hair is neg - noslipboundary_givenradiushairs.m')
        end
        
        %evaluating the Heaviside
        %hat function value
        %heaviside_fnc_value = dist_to_hair/(heaviside_width*dx);
        %sin function value
        xi = 2*dist_to_hair/(heaviside_width*dx)-1;
        heaviside_fnc_value = (1/2)*(1+xi+sin(pi*xi)/pi);
        
        regularized_fnc(outerptindex_hairs{i}(j)) = heaviside_fnc_value; 
        
    end
end
      
%setting the new velocity fields
newu = regularized_fnc.*bu;
newv = regularized_fnc.*bv;

% figure(1)    
% for i=1:num_hairs 
%     plot(xx(outerptindex_hairs{i}),yy(outerptindex_hairs{i}),'rx')
%     viscircles(hairs_center(i,:),hairs_radius+dx,'EdgeColor','g','LineWidth',1)
% end
% if (counter > 0)
%     plot(xx(outerptindex_hairs_vector(repeated_indices)),yy(outerptindex_hairs_vector(repeated_indices)),'ko')
% end

% figure(2)
% mesh(xx,yy,regularized_fnc)
% % for i=1:num_hairs
% %     viscircles(hairs_center(i,:),hairs_radius,'EdgeColor','k','LineWidth',1)
% %     viscircles(hairs_center(i,:),hairs_radius+dx,'EdgeColor','g','LineWidth',1)
% % end
% %shading(gca,'interp')

% figure(5)
% %quiver(xx,yy,bu,bv)
% skip =1;
% hold on
% quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),newu(1:skip:end,1:skip:end),newv(1:skip:end,1:skip:end))
% for i=1:num_hairs
%     viscircles(hairs_center(i,:),hairs_radius,'EdgeColor','k','LineWidth',1);
%     %viscircles(hairs_center(i,:),hairs_radius+dx,'EdgeColor','g','LineWidth',1)
% end
% 
% figure(6)
% skip =1;
% quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),bu(1:skip:end,1:skip:end),bv(1:skip:end,1:skip:end))
% hold on
% %quiver(xx,yy,newu,newv,'r')
% for i=1:num_hairs
%     viscircles(hairs_center(i,:),hairs_radius,'EdgeColor','k','LineWidth',1);
%     %viscircles(hairs_center(i,:),hairs_radius+dx,'EdgeColor','g','LineWidth',1)
% end

