function [newu,newv] = noslipboundary_notgivenradiushairs(bu,bv)

global hairs_data_filename hairs_data_filename_interior
global xshift_piv_data yshift_piv_data
global ptindex_hairs 
global x y 
 
%here we do not regularize at all but just set the velocities equal to 0
%inside the hairs that are also collecting concentration - the four
%smallest hairs are neglected
regularized_fnc = ones(size(bu)); 
[xx,yy] = ndgrid(x,y);

num_hairs = length(ptindex_hairs); 

%sets the velocity field to zero within each hair
for i=1:num_hairs
    regularized_fnc(ptindex_hairs{i}) = 0;  
end
      
%setting the new velocity fields
newu = regularized_fnc.*bu;
newv = regularized_fnc.*bv;

%allows for plotting of the hairs (based on flickdata here)
%hairs filename
load(['pivdata/' hairs_data_filename '.mat'])  
%loads the location of the hairs
hairs = eval([hairs_data_filename_interior.filename '.' hairs_data_filename_interior.hairs]); 



% figure(2)
% mesh(xx,yy,regularized_fnc)


figure(5)
skip = 1; 
%quiver(xx,yy,bu,bv)
hold on
quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),newu(1:skip:end,1:skip:end),newv(1:skip:end,1:skip:end))
for i=[1:6,8,9,11,12,14,16]
    plot(hairs_data_filename_interior.conversion_factor*hairs(i).x + xshift_piv_data(1),hairs_data_filename_interior.conversion_factor*hairs(i).y + yshift_piv_data(1),'k')
    plot(hairs_data_filename_interior.conversion_factor*hairs(i).x + xshift_piv_data(1),hairs_data_filename_interior.conversion_factor*hairs(i).y + yshift_piv_data(1),'k*')
end

figure(6)
skip = 1;
quiver(xx(1:skip:end,1:skip:end),yy(1:skip:end,1:skip:end),bu(1:skip:end,1:skip:end),bv(1:skip:end,1:skip:end))
hold on
%quiver(xx,yy,newu,newv,'r')
for i=[1:6,8,9,11,12,14,16]
     plot(hairs_data_filename_interior.conversion_factor*hairs(i).x + xshift_piv_data(1),hairs_data_filename_interior.conversion_factor*hairs(i).y + yshift_piv_data(1),'k')
     plot(hairs_data_filename_interior.conversion_factor*hairs(i).x + xshift_piv_data(1),hairs_data_filename_interior.conversion_factor*hairs(i).y + yshift_piv_data(1),'k*')
end
