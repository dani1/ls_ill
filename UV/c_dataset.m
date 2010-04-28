% This function calculates the concentrations from UV for a single dataset
function c = c_dataset ( UVpath, UVpip )

 for i = 1 : length(UVpip)
  if ~isempty(UVpip{i})

   c(i)		= UV_c(UVpath{i},UVpip{i}{:});

  end
 end

end	% c_dataset
