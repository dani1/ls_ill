% This function calculates the integral of a polynome coefficient vector
function int = polyint ( c, offset )

 d = length(c);		% polynome order

 % Change the coefficient because of integration
 for i = 1 : d
  c(end+1-i) = c(end+1-i) / i;
 end

 % Add integration factor
 if nargin == 1
  offset = 0;
 end
 int = [ c offset ];

end
