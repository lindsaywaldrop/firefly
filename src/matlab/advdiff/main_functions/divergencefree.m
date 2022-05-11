function [bu,bv] = divergencefree(xpts,ypts,bu,bv)

global piv_data_filename_interior
global uplusx_piv uminusx_piv vplusy_piv vminusy_piv
global dx dy 
global xlength ylength
global x y
        
uplusx = uplusx_piv;
uminusx = uminusx_piv;
vplusy = vplusy_piv;
vminusy = vminusy_piv;

%computing the divergence of the current velocity field
%divu = (uplusx - uminusx)/2/dx + (vplusy - vminusy)/2/dy; 
%divu2 = divergence(xpts',ypts',bu',bv');
%divu2 = divu2';
%diffdivu = divu - divu2; 
%max(max(abs(diffdivu(2:end-1,2:end-1))))
divu = divergence(xpts',ypts',bu',bv');
divu = divu';

%figure(10)
%mesh(xpts,ypts,divu)

if (piv_data_filename_interior.forcedivfree)
    
    %NOT SURE THIS WORKS

    intdivu = trapz(y,trapz(x,divu,1));
    neumannc = intdivu/2/(xlength+ylength); 

    
    NNx = size(divu);
    NNx = NNx(1);
    NNy = size(divu);
    NNy = NNy(2);
    
    %assumes dx = dy 
    [poisson_matrix_l,poisson_matrix_u,poisson_matrix] = make_poisson_matrix_noflux_noflux(divu);
    
    %setting the RHS
    divu_RHS(1:NNx*NNy,1)=0; 
    for i=1:NNx
        for j=1:NNy
            divu_RHS(i+(j-1)*NNx,1)= divu(i,j);
        end 
    end 
    %allowing neumann bc to be a constant
    for i=1:NNy 
        divu_RHS((i-1)*NNx+1)=divu_RHS((i-1)*NNx+1)-2*neumannc/dx; 
        divu_RHS(i*NNx)=divu_RHS(i*NNx)-2*neumannc/dx;;
    end
    for i=1:NNx
        divu_RHS(i) = divu_RHS(i)-2*neumannc/dy;
        divu_RHS((NNy-1)*NNx+i) = divu_RHS((NNy-1)*NNx+i)-2*neumannc/dy;
    end
    
    %choosing one point and seeting it - to deal with multiple solutions
    %for phi
    %divu_RHS(1,1) = 0;
    %divu_RHS(end,1) = 0; 
    
    phi_vector = poisson_matrix\divu_RHS;
    
    for j=1:NNy
        phi(:,j)=phi_vector((j-1)*NNx+1:j*NNx);
    end 
    
    %Note dphix = 0 on the left and right walls
    %     dphiy = 0 on the top and bottom walls 
    [dphiy,dphix] = gradient(phi,dx,dy); 
    dphix(1,:) = 0;
    dphix(end,:) = 0;
    dphiy(:,1) = 0;
    dphiy(:,end) = 0; 
    
    %divergent free velocity fields
    bu = bu - dphix;
    bv = bv - dphiy;
    
    divufinal = divergence(xpts',ypts',bu',bv');
    divufinal = divufinal';

 %   figure(11)
 %   mesh(xpts,ypts,phi)

 %   figure(12)
 %   mesh(xpts,ypts,dphix)

 %   figure(13)
 %   mesh(xpts,ypts,dphiy)

 %   figure(14)
 %   mesh(xpts,ypts,divufinal)
    
end

%--------------------------------------------------------------------------

function [matrix_l,matrix_u,poisson_matrix] = make_poisson_matrix_noflux_noflux(divu)
    
global dx dy 

%This matrix includes noflux boundary conditions in the
%x-direction and noflux in the boundary conditions in the
%y-direction

NNx = size(divu);
NNx = NNx(1);
NNy = size(divu);
NNy = NNy(2);

x_coeff=1/(dx^2);
y_coeff=1/(dy^2);

%x_coeff=1;
%y_coeff=1;

poisson_matrix = []; %sparse(NNx*NNy,NNx*NNy);

A = sparse(NNx,NNx); 

A = sparse(diag(ones(NNx,1)*-(2*x_coeff+2*y_coeff))) + ...
    sparse(diag([2;ones(NNx-2,1)]*x_coeff,1)) + ...
    sparse(diag([ones(NNx-2,1);2]*x_coeff,-1));

%full(A)

BU(1:NNx,1) = 0; 
BU(NNx+1:2*NNx,1) = 2*y_coeff; 
BU(2*NNx+1:NNx*NNy,1) = y_coeff;

BL(1:NNx*NNy-2*NNx,1) = y_coeff; 
BL(NNx*NNy-2*NNx+1:NNx*NNy-NNx,1) = 2*y_coeff;
BL(NNx*NNy-NNx+1:NNx*NNy,1) = 0; 

C_dm_diag = []; 
for i=1:NNy
  C_dm_diag = blkdiag(C_dm_diag,A); 
end

C_dm_upper = spdiags(BU,NNx,NNx*NNy,NNx*NNy);
C_dm_lower = spdiags(BL,-NNx,NNx*NNy,NNx*NNy); 

poisson_matrix = C_dm_diag + C_dm_upper + C_dm_lower;
%sets the first row to a fixed constant
%poisson_matrix(1,:) = 0;
%poisson_matrix(1,1) = 1; 
%poisson_matrix(end,:) = 0;
%poisson_matrix(end,end) = 1; 

%full(poisson_matrix)

%[matrix_l,matrix_u] = lu(poisson_matrix);
matrix_l = 0;
matrix_u = 0; 
