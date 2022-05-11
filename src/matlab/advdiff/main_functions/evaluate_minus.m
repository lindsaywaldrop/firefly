function val_minus = evaluate_minus(val, NNx, NNy, bc, xory, varargin)


if (xory == 1) %we are evaluating in the x-direction 
  i1 = 1;
  j1 = 1:NNy;
  i2 = 2:NNx;
  j2 = 1:NNy;
  i3 = NNx; 
  j3 = 1:NNy; 
  i4 = 1:NNx-1;
  j4 = 1:NNy; 
  i5 = 2;
  j5 = 1:NNy;
elseif (xory == 0) % we are working in the y-direction
  i1 = 1:NNx;
  j1 = 1; 
  i2 = 1:NNx;
  j2 = 2:NNy;
  i3 = 1:NNx;
  j3 = NNy;
  i4 = 1:NNx;
  j4 = 1:NNy-1;
  i5 = 1:NNx;
  j5 = 2;
end
  
  
  
if strcmp(bc,'periodic')
  val_minus(i1,j1)=val(i3,j3);
elseif strcmp(bc,'noflux')
  val_minus(i1,j1)=val(i5,j5);
elseif strcmp(bc,'dirichlet')
  jbd = varargin{1};  
  vel = varargin{2};
  dir_bc_value = (10*jbd)^10;
  val_minus(i1,j1) = (vel(i1,j1)<0).*dir_bc_value;  
  %val_minus(i1,j1) = dir_bc_value; 
end

val_minus(i2,j2)=val(i4,j4);
