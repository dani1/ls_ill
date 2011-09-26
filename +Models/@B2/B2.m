classdef B2
% This class describes the different models for the second virial coefficient

 properties ( Constant )

  Hard_Spheres	= inline('2/3*pi*8*R.^3*1e-27','R');		% m^3, where [R] = nm

 end

end
