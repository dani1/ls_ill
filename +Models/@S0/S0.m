classdef S0
% This class describes the different models for the structure factor at zero Q vector (which is seen y SLS)

 properties ( Constant )

  Hard_Spheres_PY	= inline('(1-phi).^4 ./ (1+2*phi).^2','phi');		% Naegele2004, p.45

 end

end
