function dcout = dcoeffvalues ( expclass, method )
% This function calculates the confint of an Experiment class
 for i = 1 : length(expclass.Point)
  tmp			= confint(expclass.Point(i).(['Fit_' method]));
  dcout(:,i)		=  0.5 * ( tmp(2,:) - tmp(1,:) );
 end
end
