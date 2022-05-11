function [velocities] = get_velocities(delt,time,explicit,flickorreturn, paths, parameters, simulation)
% 									dt_flick/2,t,explicit_vel,'flick'
%function get_velocities.m 
%input: delt, time, explicit - string stating which velocity field is desired
%output: 
%gets the velocity when desired - can include dependencies on delt and time

%global x y
%global ylength
%global handle_hairs hairs_data_filename_interior

NNx = length(simulation.x);
NNy = length(simulation.y); 

%get velocity field from functions  
velocities.u = zeros(NNx,NNy);
velocities.v = zeros(NNx,NNy); 
  
[xpts,ypts] = ndgrid(simulation.x, simulation.y); %matrix form of all the points in the bulk grid

if strcmp(explicit,'constant')
   
    %constant velocity
    bu(1:NNx,1:NNy) = -1+ypts*2/parameters.ylength; %ypts*0+1; %cos(ypts).*sin(xpts); %cos(ypts).*sin(xpts); %sin(xpts);
    bv(1:NNx,1:NNy) = ypts*0; %-sin(ypts).*cos(xpts);;   %(-(ypts-pi).^2+pi^2)/5; %sin(ypts); %-sin(ypts).*cos(xpts); %(-(ypts-pi).^2+pi^2)/5; %0; 
    
elseif strcmp(explicit, 'constant2')
    
    bu(1:NNx,1:NNy) = cos(ypts).*sin(xpts); %cos(ypts).*sin(xpts); %sin(xpts);
    bv(1:NNx,1:NNy) = -sin(ypts).*cos(xpts);;   %(-(ypts-pi).^2+pi^2)/5; %sin(ypts); %-sin(ypts).*cos(xpts); %(-(ypts-pi).^2+pi^2)/5; %0; 
    
elseif strcmp(explicit, 'constant3')
    
    bu(1:NNx,1:NNy) = 1; %cos(ypts).*sin(xpts); %cos(ypts).*sin(xpts); %sin(xpts);
    bv(1:NNx,1:NNy) = 0; %-sin(ypts).*cos(xpts);;   %(-(ypts-pi).^2+pi^2)/5; %sin(ypts); %-sin(ypts).*cos(xpts); %(-(ypts-pi).^2+pi^2)/5; %0;    
    
elseif strcmp(explicit,'constant4')
    
    if ((time+delt)<=0.5)
      bu(1:NNx,1:NNy) = -1+ypts*2/parameters.ylength;
      bv(1:NNx,1:NNy) = 0*ypts;
    elseif ((time+delt)<=1)
      bu(1:NNx,1:NNy) = (-1+ypts*2/parameters.ylength)*(2-2*(time+delt));
      bv(1:NNx,1:NNy) = 0*ypts;
    else
      bu(1:NNx,1:NNy) = 0*xpts;
      bv(1:NNx,1:NNy) = 0*ypts;
    end
    
elseif strcmp(explicit,'constant5')  
     
    if ((time+delt)<=0.5)
      bu(1:NNx,1:NNy) = -1+ypts*2/parameters.ylength;
      bv(1:NNx,1:NNy) = 0*ypts;
    elseif ((time+delt) < 1.5)
      bu(1:NNx,1:NNy) = (-1+ypts*2/parameters.ylength)*(2-2*(time+delt));
      bv(1:NNx,1:NNy) = 0*ypts;
    else
      bu(1:NNx,1:NNy) = -1*(-1+ypts*2/parameters.ylength);
      bv(1:NNx,1:NNy) = 0*ypts;
    end
    
elseif strcmp(explicit,'rigid')
    
    %Rigid Body Rotation
    bu(1:NNx,1:NNy) = -(2*pi/3)*sqrt((xpts-0.5).^2+(ypts-0.5).^2) ...
        .*sin(atan2(ypts-0.5,xpts-0.5));
    bv(1:NNx,1:NNy) = (2*pi/3)*sqrt((xpts-0.5).^2+(ypts-0.5).^2).*...
        cos(atan2(ypts-0.5,xpts-0.5));
    
elseif strcmp(explicit,'oscillating') 
    
    %stretching to ellipse and back (oscillations)
    bu(1:NNx,1:NNy) = cos(pi*(time+delt))*cos(ypts-0.5).* ...
        sin(xpts-0.5); 
    bv(1:NNx,1:NNy) = -cos(pi*(time+delt))*sin(ypts-0.5).* ...
        cos(xpts-0.5d0);
  
elseif strcmp(explicit,'rotating') 
    
    %rotating velocity
    if (mod(floor((time+delt)/0.3),2) == 0) 
      bu(1:NNx,1:NNy) = sin(pi*xpts).^2.*sin(2*pi*ypts);
      bv(1:NNx,1:NNy) = -sin(pi*ypts).^2.*sin(2*pi*xpts);
    elseif (mod(floor((time+delt)/0.3),2) == 1) 
      bu(1:NNx,1:NNy) = -sin(pi*xpts).^2.*sin(2*pi*ypts);
      bv(1:NNx,1:NNy) = sin(pi*ypts).^2.*sin(2*pi*xpts);
    else
      fprintf('velocity error')
    end 
    
elseif strcmp(explicit,'piv_data')
    
        %read in data
        [velocities] = read_in_velocity_data_p2(xpts, ypts, flickorreturn, paths, parameters);   
    
%         if (handle_hairs)
%             if (hairs_data_filename_interior.givenradius)
%                 %add in no slip boundary condition if hairs are small 
%                 [bu,bv] = noslipboundary_givenradiushairs(bu,bv); 
%             elseif (~hairs_data_filename_interior.givenradius)
%                 [bu,bv] = noslipboundary_notgivenradiushairs(bu,bv); 
%             end
%             appending_vel_piv_bcs(bu,bv); 
%         end
    
        %[bu,bv] = divergencefree(xpts,ypts,bu,bv); 
        
else 
    error('explicit velocity is not a valid choice!'); 
  
end

