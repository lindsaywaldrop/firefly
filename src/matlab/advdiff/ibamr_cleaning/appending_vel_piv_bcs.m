function appending_vel_piv_bcs(bu,bv) 
global uplusx_piv uminusx_piv uplusy_piv uminusy_piv
global vplusx_piv vminusx_piv vplusy_piv vminusy_piv
global uplus2x_piv uminus2x_piv vplus2y_piv vminus2y_piv

%since we modified the interior of the vel fields to deal with the hairs
%this must also be done for these velocities

uplusx_piv(1:end-1,:) = bu(2:end,:);
uplus2x_piv(1:end-1,:) = uplusx_piv(2:end,:);
uminusx_piv(2:end,:) = bu(1:end-1,:);
uminus2x_piv(2:end,:) = uminusx_piv(1:end-1,:);
uplusy_piv(:,1:end-1) = bu(:,2:end);
uminusy_piv(:,2:end) = bu(:,1:end-1); 


vplusx_piv(1:end-1,:) = bv(2:end,:);
vminusx_piv(2:end,:) = bv(1:end-1,:);
vplusy_piv(:,1:end-1) = bv(:,2:end);
vplus2y_piv(:,1:end-1) = vplusy_piv(:,2:end);
vminusy_piv(:,2:end) = bv(:,1:end-1); 
vminus2y_piv(:,2:end) = vminusy_piv(:,1:end-1); 