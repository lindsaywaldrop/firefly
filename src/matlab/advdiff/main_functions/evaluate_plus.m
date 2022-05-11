function val_plus = evaluate_plus(val, NNx, NNy, bc, xory, varargin)

if (xory==1) % we are working in the x-direction
  i1 = 1;
  j1 = 1:NNy;
  i2 = 2:NNx;
  j2 = 1:NNy;
  i3 = NNx; 
  j3 = 1:NNy; 
  i4 = 1:NNx-1;
  j4 = 1:NNy; 
  i5 = NNx-1;
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
  j5 = NNy-1;
end

if strcmp(bc,'periodic')
  val_plus(i3,j3) = val(i1,j1); 
elseif strcmp(bc,'noflux')
  val_plus(i3,j3) = val(i5,j5); 
elseif strcmp(bc,'dirichlet')
  jbd = varargin{1};  
  vel = varargin{2};
  if (length(varargin) == 3)
      cdbc = varargin{3}; 
  else
      cdbc = 0; 
  end
  dir_bc_value = (10*jbd)^10;
  %val_plus(i3,j3) = (vel(i3,j3)>0).*dir_bc_value; 
  if (vel(i3,j3)>0)
      val_plus(i3,j3) = dir_bc_value; 
  else 
     
     val_plus(i3,j3) = cdbc; 
  end
  %val_plus(i3,j3) = dir_bc_value; 
end
  

val_plus(i4,j4) = val(i2,j2); 


  

