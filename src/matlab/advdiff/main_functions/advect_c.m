function [simulation] = advect_c(delt, bc, method, parameters, simulation, velocities)

%global dx dy c u v
%global uplusx_piv uminusx_piv uplusy_piv uminusy_piv
%global vplusx_piv vminusx_piv vplusy_piv vminusy_piv
%global uplus2x_piv uminus2x_piv vplus2y_piv vminus2y_piv
%global cplusx_dbc

NNx = size(simulation.c);
NNx = NNx(1);
NNy = size(simulation.c);
NNy = NNy(2); 

%boundary conditions
%note for WENO we also reset bcs
if strcmp(bc,'dirichlet')
    %Dirichlet boundary conditions for c
    if strcmp(method,'upwind')
        
        uplusx = velocities.uplusx_piv;
        uminusx = velocities.uminusx_piv;
        vplusy = velocities.vplusy_piv;
        vminusy = velocities.vminusy_piv;
        
        cplusx = evaluate_plus(simulation.c, NNx, NNy, 'dirichlet', 1, 1, velocities.u, parameters.cplusx_dbc);
        cminusx = evaluate_minus(simulation.c, NNx, NNy, 'dirichlet', 1, 1, velocities.u);
        cplusy = evaluate_plus(simulation.c, NNx, NNy, 'dirichlet', 0, 1, velocities.v);
        cminusy = evaluate_minus(simulation.c, NNx, NNy, 'dirichlet', 0, 1, velocities.v);
        
    elseif strcmp(method,'LaxWendroff')
        %NEEDS FIXING
        %needed for Lax-Wendroff
        %cplusxplusy
        cplusxplusy(1:NNx-1, 1:NNy-1) = simulation.c(2:NNx,2:NNy);
        cplusxplusy(NNx, 1:NNy-1) = 0;
        cplusxplusy(NNx, NNy) = 0;
        cplusxplusy(1:NNx-1, NNy) = 0;
        %cplusxminusy
        cplusxminusy(1:NNx-1, 2:NNy) = simulation.c(2:NNx,1:NNy-1);
        cplusxminusy(NNx, 2:NNy) = 0;
        cplusxminusy(NNx, 1) = 0;
        cplusxminusy(1:NNx-1, 1) = 0;
        %cminusxplusy
        cminusxplusy(2:NNx, 1:NNy-1) = simulation.c(1:NNx-1,2:NNy);
        cminusxplusy(1, 1:NNy-1) = 0;
        cminusxplusy(1, NNy) = 0;
        cminusxplusy(2:NNx, NNy) = 0;
        %cminusxminusy
        cminusxminusy(2:NNx, 2:NNy) = simulation.c(1:NNx-1,1:NNy-1);
        cminusxminusy(1, 2:NNy) = 0;
        cminusxminusy(1, 1) = 0;
        cminusxminusy(2:NNx, 1) = 0;

        %boundary conditions on u and v using extrapolation based on piv data
        uplusx = velocities.uplusx_piv;
        uminusx = velocities.uminusx_piv;
        vplusx = velocities.vplusx_piv;
        vminusx = velocities.vminusx_piv;
        uplusy = velocities.uplusy_piv;
        uminusy = velocities.uminusy_piv;
        vplusy = velocities.vplusy_piv;
        vminusy = velocities.vminusy_piv;
        
    elseif strcmp(method,'weno')
        %needed for WENO 
        uplusx = velocities.uplusx_piv;
        uminusx = velocities.uminusx_piv;
        uplus2x = velocities.uplus2x_piv;
        uminus2x = velocities.uminus2x_piv;
        vplusy = velocities.vplusy_piv;
        vplus2y = velocities.vplus2y_piv;
        vminusy = velocities.vminusy_piv;
        vminus2y = velocities.vminus2y_piv;
        
        cplusx = evaluate_plus(simulation.c, NNx, NNy, 'dirichlet', 1, 1, velocities.uplusx_piv, parameters.cplusx_dbc);
        cminusx = evaluate_minus(simulation.c, NNx, NNy, 'dirichlet', 1, 1, velocities.uminusx_piv);
        cplusy = evaluate_plus(simulation.c, NNx, NNy, 'dirichlet', 0, 1, velocities.vplusy_piv);
        cminusy = evaluate_minus(simulation.c, NNx, NNy, 'dirichlet', 0, 1, velocities.vminusy_piv);
        
        cplus2x = evaluate_plus(cplusx, NNx, NNy, 'dirichlet', 1, 2, velocities.uplus2x_piv, parameters.cplusx_dbc);
        cminus2x = evaluate_minus(cminusx, NNx, NNy, 'dirichlet', 1, 2, velocities.uminus2x_piv);
        cplus2y = evaluate_plus(cplusy, NNx, NNy, 'dirichlet', 0, 2, velocities.vplus2y_piv);
        cminus2y = evaluate_minus(cminusy, NNx, NNy, 'dirichlet', 0, 2, velocities.vminus2y_piv);
               
    end
    
elseif strcmp(bc,'periodic_noflux')
    %Previously with no flux and periodic boundary conditions 
    %Eventually want to add this in as a function! 
    %is a little weird for WENO - see below
    
    %unext = u;
    %vnext = v; 
    
    cplusx = evaluate_plus(simulation.c, NNx, NNy, 'periodic', 1); 
    cminusx = evaluate_minus(simulation.c, NNx, NNy, 'periodic', 1); 
    cplusy = evaluate_plus(simulation.c, NNx, NNy, 'noflux', 0);
    cminusy = evaluate_minus(simulation.c, NNx, NNy, 'noflux', 0);

    if strcmp(method,'LaxWendroff')
        %needed for Lax-Wendroff    
        cplusxplusy(1:NNx-1, 1:NNy-1) = simulation.c(2:NNx, 2:NNy);
        cplusxplusy(NNx, 1:NNy-1) =simulation.c(1, 2:NNy);
        cplusxplusy(NNx, NNy) = simulation.c(1, NNy-1);
        cplusxplusy(1:NNx-1, NNy) = simulation.c(2:NNx, NNy-1);

        cplusxminusy(1:NNx-1, 2:NNy) = simulation.c(2:NNx, 1:NNy-1);
        cplusxminusy(NNx, 2:NNy) = simulation.c(1, 1:NNy-1);
        cplusxminusy(NNx, 1) = simulation.c(1, 2);
        cplusxminusy(1:NNx-1, 1) = simulation.c(2:NNx, 2);

        cminusxplusy(2:NNx, 1:NNy-1) = simulation.c(1:NNx-1, 2:NNy);
        cminusxplusy(1, 1:NNy-1) = simulation.c(NNx, 2:NNy);
        cminusxplusy(1, NNy) = simulation.c(NNx, NNy-1);
        cminusxplusy(2:NNx, NNy) = simulation.c(1:NNx-1, NNy-1);

        cminusxminusy(2:NNx, 2:NNy) = simulation.c(1:NNx-1, 1:NNy-1);
        cminusxminusy(1, 2:NNy) = simulation.c(NNx, 1:NNy-1);
        cminusxminusy(1, 1) = simulation.c(NNx, 2);
        cminusxminusy(2:NNx, 1) = simulation.c(1:NNx-1, 2);
 
        uplusx = evaluate_plus(velocities.u, NNx, NNy, 'periodic', 1);
        uminusx = evaluate_minus(velocities.u, NNx, NNy, 'periodic', 1);
        vplusx = evaluate_plus(velocities.v, NNx, NNy, 'periodic', 1);
        vminusx = evaluate_minus(velocities.v, NNx, NNy, 'periodic', 1);

        %IS THIS CORRECT? 
        vplusy(:, 1:NNy-1) = velocities.v(:, 2:NNy);
        vplusy(:, NNy) = -velocities.v(:, NNy-1);                
        vminusy(:, 2:NNy)=v(:, 1:NNy-1);
        vminusy(:, 1)= -velocities.v(:, 2);

        uplusy = evaluate_plus(velocities.u, NNx, NNy, 'noflux', 0);
        uminusy = evaluate_minus(velocities.u, NNx, NNy, 'noflux', 0);  
        
    elseif strcmp(method,'weno')
        %needed for WENO 
        
        %boundary conditions on u and v using extrapolation based on piv
        %data - should be something different if not constant velocity
        uplusx = velocities.uplusx_piv;
        uminusx = velocities.uminusx_piv;
        uplus2x = velocities.uplus2x_piv;
        uminus2x = velocities.uminus2x_piv;
        vplusx = velocities.vplusx_piv;
        vminusx = velocities.vminusx_piv;
        uplusy = velocities.uplusy_piv;
        uminusy = velocities.uminusy_piv;
        vplusy = velocities.vplusy_piv;
        vminusy = velocities.vminusy_piv;
        vplus2y = velocities.vplus2y_piv;
        vminus2y = velocities.vminus2y_piv;
    
        cplus2x = evaluate_plus(cplusx, NNx, NNy, 'periodic', 1); 
        cminus2x = evaluate_minus(cminusx, NNx, NNy, 'periodic', 1); 
        %not correct if c is not constant in y direction 
        cplus2y = evaluate_plus(cplusy, NNx, NNy, 'noflux', 0);
        cminus2y = evaluate_minus(cminusy, NNx, NNy, 'noflux', 0);
               
    end
    
else
    error('no such bc in advect_c.m')
end

if strcmp(method,'LaxWendroff')
    %using a second order method (Lax-Wendroff) to update c
    %this does not work with Dirichlet 0 boundary conditions because it causes
    %buildup near the wall - maybe do some kind of extrapolation of c
    %instead...
    %DOES NOT SOLVE THE CONSERVATIVE FORM OF THE EQUATIONS
    c_advx = (-delt*velocities.u.*(cplusx-cminusx)/2/parameters.dx)...
        +(delt^2*velocities.u.*(((uplusx-uminusx)/2/parameters.dx).*((cplusx-cminusx)/2/parameters.dx))/2)...
        +(delt^2*velocities.u.*(velocities.u.*(cplusx-2*c+cminusx)/parameters.dx^2)/2);
    
    c_advy = (-delt*velocities.v.*(cplusy-cminusy)/2/parameters.dy)...
        +(delt^2*velocities.v.*(((vplusy-vminusy)/2/parameters.dy).*((cplusy-cminusy)/2/parameters.dy))/2)...
        +(delt^2*velocities.v.*(velocities.v.*(cplusy-2*c+cminusy)/parameters.dy^2)/2);
    
    c_advxy = (delt^2*velocities.u.*(((vplusx-vminusx)/2/parameters.dx).*((cplusy-cminusy)/2/parameters.dy))/2)...
        +(delt^2*velocities.v.*(((uplusy-uminusy)/2/parameters.dy).*((cplusx-cminusx)/2/parameters.dx))/2)...
        +(delt^2*velocities.u.*velocities.v.*(cplusxplusy-cplusxminusy-cminusxplusy+cminusxminusy)/4/parameters.dx/parameters.dy);
    
    simulation.c = simulation.c + c_advx + c_advy + c_advxy;
    
    %former version
    % %using a second order method (Lax-Wendroff)                                               
    % c_advx=(-delt*u.*(cplusx-cminusx)/2/dx)...
    %        +(delt^2*u.*((((uplusx-uminusx)/2/dx).*((cplusx-cminusx)/2/dx))...
    %                         +(u.*(cplusx-2*c+cminusx)/dx^2))/2)...
    %        -(delt*(unext-u).*(cplusx-cminusx)/4/dx);
    % 
    % c_advy=(-delt*v.*(cplusy-cminusy)/2/dy)...
    %        +(delt^2*v.*((((vplusy-vminusy)/2/dy).*((cplusy-cminusy)/2/dy))...
    %                         +(v.*(cplusy-2*c+cminusy)/dy^2))/2)...
    %        -(delt*(vnext-v).*(cplusy-cminusy)/4/dy);
    % 
    % c_advxy=delt^2*((u.*(vplusx-vminusx).*(cplusy-cminusy))+...
    %                 (v.*(uplusy-uminusy).*(cplusx-cminusx)))/8/dx/dy+...
    %         delt^2*u.*v.*(cplusxplusy-cplusxminusy- ...
    %                               cminusxplusy+cminusxminusy)/4/dx/dy; 
    % c=c+c_advx+c_advy+c_advxy;
elseif strcmp(method,'upwind')
    %using a first order method (Upwind) to update c 
    %c_advx = -delt*max(u,0).*(c-cminusx)/dx-delt*min(u,0).*(cplusx-c)/dx;
    %c_advy = -delt*max(v,0).*(c-cminusy)/dy-delt*min(v,0).*(cplusy-c)/dy;

    c_advx = -delt*(velocities.u>0).*(velocities.u.*simulation.c-uminusx.*cminusx)/...
    		parameters.dx-delt*(velocities.u<0).*(uplusx.*cplusx-velocities.u.*simulation.c)/parameters.dx;
    c_advy = -delt*(velocities.v>0).*(velocities.v.*simulation.c-vminusy.*cminusy)/...
    		parameters.dy-delt*(velocities.v<0).*(vplusy.*cplusy-velocities.v.*simulation.c)/parameters.dy;
    
    simulation.c = simulation.c + c_advx+c_advy;
    
elseif strcmp(method,'weno')
    %WENO with k=2 from Shu 1997
    
    %timestepping - TVD RK3 from Shu 1997
    
    c1 = simulation.c + delt*wenoflux(parameters, simulation.c, velocities.u, velocities.v,...
    									uplusx,uminusx,uplus2x,uminus2x,...
    									vplusy,vminusy,vplus2y,vminus2y,cplusx,...
    									cminusx,cplusy,cminusy,cplus2x,cminus2x,...
    									cplus2y,cminus2y);
                                    
    %setting boundary conditions for c1
    %DIRICHLET
    cplusx = evaluate_plus(c1, NNx, NNy, 'dirichlet', 1, 1, uplusx, parameters.cplusx_dbc);
    cminusx = evaluate_minus(c1, NNx, NNy, 'dirichlet', 1, 1, uminusx);
    cplusy = evaluate_plus(c1, NNx, NNy, 'dirichlet', 0, 1, vplusy);
    cminusy = evaluate_minus(c1, NNx, NNy, 'dirichlet', 0, 1, vminusy);
    cplus2x = evaluate_plus(cplusx, NNx, NNy, 'dirichlet', 1, 2, uplus2x, parameters.cplusx_dbc);
    cminus2x = evaluate_minus(cminusx, NNx, NNy, 'dirichlet', 1, 2, uminus2x);
    cplus2y = evaluate_plus(cplusy, NNx, NNy, 'dirichlet', 0, 2, vplus2y);
    cminus2y = evaluate_minus(cminusy, NNx, NNy, 'dirichlet', 0, 2, vminus2y);
    
    %not correct if c is not constant in y direction
    %PERIODIC/NOFLUX
    %cplusx = evaluate_plus(c1,NNx,NNy,'periodic',1); 
    %cminusx = evaluate_minus(c1,NNx,NNy,'periodic',1); 
    %cplusy = evaluate_plus(c1,NNx,NNy,'noflux',0);
    %cminusy = evaluate_minus(c1,NNx,NNy,'noflux',0);
    %cplus2x = evaluate_plus(cplusx,NNx,NNy,'periodic',1); 
    %cminus2x = evaluate_minus(cminusx,NNx,NNy,'periodic',1); 
    %cplus2y = evaluate_plus(cplusy,NNx,NNy,'noflux',0);
    %cminus2y = evaluate_minus(cminusy,NNx,NNy,'noflux',0);
    
    c2 = 3*simulation.c/4 + c1/4 + delt*wenoflux(parameters, c1, velocities.u, velocities.v, uplusx,...
    												uminusx, uplus2x, uminus2x, vplusy, vminusy,... 
    												vplus2y, vminus2y, cplusx, cminusx, cplusy,... 
    												cminusy, cplus2x, cminus2x, cplus2y, cminus2y)/4;
    
    %setting boundary conditions for c2
    %DIRICHLET
    cplusx = evaluate_plus(c2, NNx, NNy, 'dirichlet', 1, 1, uplusx, parameters.cplusx_dbc);
    cminusx = evaluate_minus(c2, NNx, NNy, 'dirichlet', 1, 1, uminusx);
    cplusy = evaluate_plus(c2, NNx, NNy, 'dirichlet', 0, 1, vplusy);
    cminusy = evaluate_minus(c2, NNx, NNy, 'dirichlet', 0, 1, vminusy);
    cplus2x = evaluate_plus(cplusx, NNx, NNy, 'dirichlet', 1, 2, uplus2x, parameters.cplusx_dbc);
    cminus2x = evaluate_minus(cminusx, NNx, NNy, 'dirichlet', 1, 2, uminus2x);
    cplus2y = evaluate_plus(cplusy, NNx, NNy, 'dirichlet', 0, 2, vplus2y);
    cminus2y = evaluate_minus(cminusy, NNx, NNy, 'dirichlet', 0, 2, vminus2y);

    %not correct if c is not constant in y direction
    %PERIODIC/NOFLUX
    %cplusx = evaluate_plus(c2,NNx,NNy,'periodic',1); 
    %cminusx = evaluate_minus(c2,NNx,NNy,'periodic',1); 
    %cplusy = evaluate_plus(c2,NNx,NNy,'noflux',0);
    %cminusy = evaluate_minus(c2,NNx,NNy,'noflux',0);
    %cplus2x = evaluate_plus(cplusx,NNx,NNy,'periodic',1); 
    %cminus2x = evaluate_minus(cminusx,NNx,NNy,'periodic',1); 
    %cplus2y = evaluate_plus(cplusy,NNx,NNy,'noflux',0);
    %cminus2y = evaluate_minus(cminusy,NNx,NNy,'noflux',0);
    
    simulation.c = simulation.c/3 + 2*c2/3 + 2*delt*...
    				wenoflux(parameters, c2, velocities.u, velocities.v, uplusx, uminusx, uplus2x,...
    				uminus2x, vplusy, vminusy, vplus2y, vminus2y, cplusx, cminusx,...
    				cplusy, cminusy, cplus2x, cminus2x, cplus2y, cminus2y)/3;
  

    %c = c + delt*wenoflux(c,u,v,uplusx,uminusx,uplus2x,uminus2x,vplusy,vminusy,vplus2y,vminus2y,cplusx,cminusx,cplusy,cminusy,cplus2x,cminus2x,cplus2y,cminus2y);
    
else
   error('no such method in advect_c.m') 
end


function [f] = wenoflux(parameters, cc, uu, vv, uuplusx, uuminusx, uuplus2x, uuminus2x, ...
						vvplusy, vvminusy, vvplus2y, vvminus2y, ccplusx, ccminusx, ccplusy,...
						ccminusy, ccplus2x, ccminus2x, ccplus2y, ccminus2y)
%takes care of spatial derivatives for weno

%Note - uu,vv,uuplus,uuminus,vvplus,and vvminus should not be changing
%over time

%global dx parameters.dy weno_eps

NNx = size(cc);
NNx = NNx(1);
NNy = size(cc);
NNy = NNy(2); 

%r = 0
fxhalf0minus(1,:) = uuminusx(1,:).*ccminusx(1,:)/2 + uu(1,:).*cc(1,:)/2;
fxhalf0minus(2:NNx+1,:) = uu.*cc/2 + uuplusx.*ccplusx/2;
fxhalf0plus(1,:) = 3*uu(1,:).*cc(1,:)/2 - uuplusx(1,:).*ccplusx(1,:)/2;
fxhalf0plus(2:NNx+1,:) = 3*uuplusx.*ccplusx/2 - uuplus2x.*ccplus2x/2;

fyhalf0minus(:,1) = vvminusy(:,1).*ccminusy(:,1)/2 + vv(:,1).*cc(:,1)/2;
fyhalf0minus(:,2:NNy+1) = vv.*cc/2 + vvplusy.*ccplusy/2;
fyhalf0plus(:,1) = 3*vv(:,1).*cc(:,1)/2 - vvplusy(:,1).*ccplusy(:,1)/2;
fyhalf0plus(:,2:NNy+1) = 3*vvplusy.*ccplusy/2 - vvplus2y.*ccplus2y/2;

%r=1
fxhalf1minus(1,:) = -uuminus2x(1,:).*ccminus2x(1,:)/2 + 3*uuminusx(1,:).*ccminusx(1,:)/2;
fxhalf1minus(2:NNx+1,:) = -uuminusx.*ccminusx/2 + 3*uu.*cc/2;
fxhalf1plus(1,:) = uuminusx(1,:).*ccminusx(1,:)/2 + uu(1,:).*cc(1,:)/2;
fxhalf1plus(2:NNx+1,:) = uu.*cc/2 + uuplusx.*ccplusx/2;

fyhalf1minus(:,1) = -vvminus2y(:,1).*ccminus2y(:,1)/2 + 3*vvminusy(:,1).*ccminusy(:,1)/2;
fyhalf1minus(:,2:NNy+1) = -vvminusy.*ccminusy/2 + 3*vv.*cc/2;
fyhalf1plus(:,1) = vvminusy(:,1).*ccminusy(:,1)/2 + vv(:,1).*cc(:,1)/2;
fyhalf1plus(:,2:NNy+1) = vv.*cc/2 + vvplusy.*ccplusy/2;

%weights
epsilon = parameters.weno_eps; 

d0 = 2/3; 
d1 = 1/3;

%note since we do not sqaure the betas here -> we do in the alphas
betax0minus(1,:) = uu(1,:).*cc(1,:) - uuminusx(1,:).*ccminusx(1,:);
betax0minus(2:NNx+1,:) = uuplusx.*ccplusx - uu.*cc;
betax1minus(1,:) = uuminusx(1,:).*ccminusx(1,:) - uuminus2x(1,:).*ccminus2x(1,:);
betax1minus(2:NNx+1,:) = uu.*cc - uuminusx.*ccminusx; 
betax0plus(1,:) = uuplusx(1,:).*ccplusx(1,:) - uu(1,:).*cc(1,:);
betax0plus(2:NNx+1,:) = uuplus2x.*ccplus2x - uuplusx.*ccplusx;
betax1plus = betax0minus; 

betay0minus(:,1) = vv(:,1).*cc(:,1) - vvminusy(:,1).*ccminusy(:,1);
betay0minus(:,2:NNy+1) = vvplusy.*ccplusy - vv.*cc;
betay1minus(:,1) = vvminusy(:,1).*ccminusy(:,1) - vvminus2y(:,1).*ccminus2y(:,1);
betay1minus(:,2:NNy+1) = vv.*cc - vvminusy.*ccminusy; 
betay0plus(:,1) = vvplusy(:,1).*ccplusy(:,1) - vv(:,1).*cc(:,1);
betay0plus(:,2:NNy+1) = vvplus2y.*ccplus2y - vvplusy.*ccplusy;
betay1plus = betay0minus; 

alphax0minus = d0./((epsilon + betax0minus.^2).^2);
alphax1minus = d1./((epsilon + betax1minus.^2).^2);
alphax0plus = d1./((epsilon + betax0plus.^2).^2);
alphax1plus = d0./((epsilon + betax1plus.^2).^2);

alphay0minus = d0./((epsilon + betay0minus.^2).^2);
alphay1minus = d1./((epsilon + betay1minus.^2).^2);
alphay0plus = d1./((epsilon + betay0plus.^2).^2);
alphay1plus = d0./((epsilon + betay1plus.^2).^2);

omegax0minus = alphax0minus./(alphax0minus + alphax1minus);
omegax1minus = alphax1minus./(alphax0minus + alphax1minus);
omegax0plus = alphax0plus./(alphax0plus + alphax1plus);
omegax1plus = alphax1plus./(alphax0plus + alphax1plus);

omegay0minus = alphay0minus./(alphay0minus + alphay1minus);
omegay1minus = alphay1minus./(alphay0minus + alphay1minus);
omegay0plus = alphay0plus./(alphay0plus + alphay1plus);
omegay1plus = alphay1plus./(alphay0plus + alphay1plus);


%putting it all together
fxhalfminus = omegax0minus.*fxhalf0minus + omegax1minus.*fxhalf1minus;
fxhalfplus = omegax0plus.*fxhalf0plus + omegax1plus.*fxhalf1plus;

fyhalfminus = omegay0minus.*fyhalf0minus + omegay1minus.*fyhalf1minus;
fyhalfplus = omegay0plus.*fyhalf0plus + omegay1plus.*fyhalf1plus;

ahalfx(1,:) = betax0minus(1,:)./(cc(1,:)-ccminusx(1,:));
ahalfx(2:NNx+1,:) = betax0minus(2:end,:)./(ccplusx-cc);
ahalfy(:,1) = betay0minus(:,1)./(cc(:,1)-ccminusy(:,1));
ahalfy(:,2:NNy+1) = betay0minus(:,2:end)./(ccplusy-cc);

fxhalf = (ahalfx >= 0).*fxhalfminus + (ahalfx < 0).*fxhalfplus;
fyhalf =  (ahalfy >= 0).*fyhalfminus + (ahalfy < 0).*fyhalfplus;

fx = (fxhalf(1:end-1,:) - fxhalf(2:end,:))/parameters.dx;
fy = (fyhalf(:,1:end-1) - fyhalf(:,2:end))/parameters.dy;

f = fx+fy;  


