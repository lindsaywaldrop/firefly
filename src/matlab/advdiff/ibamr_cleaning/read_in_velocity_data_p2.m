function [velocities] = read_in_velocity_data_p2(xpts, ypts, flickorreturn, paths, parameters)
%function read_in_velocity_data_p2.m 
%input: xpts, ypts - grid for velocity field
%output: bu, bv - velocity field on xpts and ypts
%sets u and v based on the experimental piv data 
%interpolating to the c grid

%global pathbase_piv
%global piv_data_filename piv_data_returnfilename piv_data_filename_interior 
%global xshift_piv_data yshift_piv_data

%global uplusx_piv uminusx_piv uplusy_piv uminusy_piv
%global vplusx_piv vminusx_piv vplusy_piv vminusy_piv
%global uplus2x_piv uminus2x_piv vplus2y_piv vminus2y_piv
%global dx dy

if strcmp(flickorreturn,'rest')
    velocities.u = xpts*0;
    velocities.v = ypts*0;  
    
    uplusx_piv = velocities.u;
    uplus2x_piv = velocities.u;
    uminusx_piv = velocities.u;
    uminus2x_piv = velocities.u;
    uplusy_piv = velocities.u;
    uminusy_piv = velocities.u;


    vplusx_piv = velocities.v;
    vminusx_piv = velocities.v;
    vplusy_piv = velocities.v;
    vplus2y_piv = velocities.v;
    vminusy_piv = velocities.v;
    vminus2y_piv = velocities.v;
    
else

    %filename
    if strcmp(flickorreturn,'flick') 
        flickdata = load(strcat(paths.pathbase_piv, parameters.piv_data_filename, '.mat'));
        shift_comp = 1; 
    %elseif strcmp(flickorreturn,'return')
    %    load(['pivdata/', piv_data_returnfilename '.mat']) 
    %    shift_comp = 2; 
    else
        error('not a type of piv velocity')
    end
    
    %data is in cm and needs to be transposed - NO LONGER TRUE 
    xpiv = parameters.piv_data_filename_interior.conversion_factor*(eval(['flickdata.' parameters.piv_data_filename_interior.x]));      
    ypiv = parameters.piv_data_filename_interior.conversion_factor*(eval(['flickdata.' parameters.piv_data_filename_interior.y]));

    %shift the data so that it agrees with the concentration grid
    %minus 2 allows for extrap for weno2
    xpiv = xpiv + parameters.xshift_piv_data(shift_comp);
    ypiv = ypiv + parameters.yshift_piv_data(shift_comp); 

    %if data is in different units and needs to be transposed
    upiv = parameters.piv_data_filename_interior.conversion_factor*(eval(['flickdata.' parameters.piv_data_filename_interior.u]));      
    vpiv = parameters.piv_data_filename_interior.conversion_factor*(eval(['flickdata.' parameters.piv_data_filename_interior.v]));

    %temp if there are nans in the original data
    [iupiv] = find(isnan(upiv)==1);
    upiv(iupiv) = 0; 
    [ivpiv] = find(isnan(vpiv)==1);
    vpiv(ivpiv) = 0;

    %interpolating experimental data to the c grid 
    velocities.u = griddata(xpiv',ypiv',upiv',xpts,ypts);
    velocities.v = griddata(xpiv',ypiv',vpiv',xpts,ypts);

    %extrapolates u and v experimental data for when it is needed in advect_c

    %extrapolation based on values at boundary 
    % uplusx_piv_vec = bu(end,:);
    % uminusx_piv_vec = bu(1,:); 
    % uplus2x_piv_vec = bu(end,:);
    % uminus2x_piv_vec = bu(1,:); 
    % uplusy_piv_vec = bu(:,end);
    % uminusy_piv_vec = bu(:,1); 
    % vplusx_piv_vec = bv(end,:); 
    % vminusx_piv_vec = bv(1,:); 
    % vplusy_piv_vec = bv(:,end);
    % vminusy_piv_vec = bv(:,1); 
    % vplus2y_piv_vec = bv(:,end);
    % vminus2y_piv_vec = bv(:,1); 

    %extrapolation (linear)
	uplusx_piv_vec = griddata(xpiv',ypiv',upiv',xpts(end,:)+parameters.dx,ypts(end,:));
    uminusx_piv_vec = griddata(xpiv',ypiv',upiv',xpts(1,:)-parameters.dx,ypts(1,:));
    uplus2x_piv_vec = griddata(xpiv',ypiv',upiv',xpts(end,:)+2*parameters.dx,ypts(end,:));
    uminus2x_piv_vec = griddata(xpiv',ypiv',upiv',xpts(1,:)-2*parameters.dx,ypts(1,:));
    uplusy_piv_vec = griddata(xpiv',ypiv',upiv',xpts(:,end),ypts(:,end)+parameters.dy);
    uminusy_piv_vec = griddata(xpiv',ypiv',upiv',xpts(:,1),ypts(:,1)-parameters.dy);
    vplusx_piv_vec = griddata(xpiv',ypiv',vpiv',xpts(end,:)+parameters.dx,ypts(end,:));
    vminusx_piv_vec = griddata(xpiv',ypiv',vpiv',xpts(1,:)-parameters.dx,ypts(1,:));
    vplusy_piv_vec = griddata(xpiv',ypiv',vpiv',xpts(:,end),ypts(:,end)+parameters.dy);
    vminusy_piv_vec = griddata(xpiv',ypiv',vpiv',xpts(:,1),ypts(:,1)-parameters.dy);
    vplus2y_piv_vec = griddata(xpiv',ypiv',vpiv',xpts(:,end),ypts(:,end)+2*parameters.dy);
    vminus2y_piv_vec = griddata(xpiv',ypiv',vpiv',xpts(:,1),ypts(:,1)-2*parameters.dy);

    %different method of interpolation
    %uplusx_piv_vec = interp2(xpts',ypts',bu',xpts(end,:)+parameters.dx,ypts(end,:),'spline');
    %uminusx_piv_vec = interp2(xpts',ypts',bu',xpts(1,:)-parameters.dx,ypts(1,:),'spline');
    %uplusy_piv_vec = interp2(xpts',ypts',bu',xpts(:,end),ypts(:,end)+parameters.dy,'spline');
    %uminusy_piv_vec = interp2(xpts',ypts',bu',xpts(:,1),ypts(:,1)-parameters.dy,'spline');
    %vplusx_piv_vec = interp2(xpts',ypts',bv',xpts(end,:)+parameters.dx,ypts(end,:),'spline');
    %vminusx_piv_vec = interp2(xpts',ypts',bv',xpts(1,:)-parameters.dx,ypts(1,:),'spline');
    %vplusy_piv_vec = interp2(xpts',ypts',bv',xpts(:,end),ypts(:,end)+parameters.dy,'spline');
    %vminusy_piv_vec = interp2(xpts',ypts',bv',xpts(:,1),ypts(:,1)-parameters.dy,'spline');


    velocities.uplusx_piv = [velocities.u(2:end,:);uplusx_piv_vec];
    velocities.uplus2x_piv = [velocities.uplusx_piv(2:end,:);uplus2x_piv_vec];
    velocities.uminusx_piv = [uminusx_piv_vec;velocities.u(1:end-1,:)];
    velocities.uminus2x_piv = [uminus2x_piv_vec;velocities.uminusx_piv(1:end-1,:)];
    velocities.uplusy_piv = [velocities.u(:,2:end) uplusy_piv_vec];
    velocities.uminusy_piv = [uminusy_piv_vec velocities.u(:,1:end-1)]; 


    velocities.vplusx_piv = [velocities.v(2:end,:);vplusx_piv_vec];
    velocities.vminusx_piv = [vminusx_piv_vec;velocities.v(1:end-1,:)];
    velocities.vplusy_piv = [velocities.v(:,2:end) vplusy_piv_vec];
    velocities.vplus2y_piv = [velocities.vplusy_piv(:,2:end) vplus2y_piv_vec];
    velocities.vminusy_piv = [vminusy_piv_vec velocities.v(:,1:end-1)]; 
    velocities.vminus2y_piv = [vminus2y_piv_vec velocities.vminusy_piv(:,1:end-1)]; 
    
end
