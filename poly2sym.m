% Poly2sym function, written by foot
function pol = poly2sym( c )
% c are the polynome coefficients
 
 d = length(c);		% polynome order

 pol = '';

 for i = 1 : d
  if c(end+1-i) >= 0
   s = '+';
  else
   s = '-';
  end

  if i == 1
   pol = [s ' ' num2str(abs(c(end+1-i)))];
  elseif i == 2
   pol = [s ' ' num2str(abs(c(end+1-i))) ' * x ' pol];
  else
   pol = [s ' ' num2str(abs(c(end+1-i))) ' * x.^' num2str(i-1) ' ' pol];
  end

 end

 pol = inline(pol,'x');

end
