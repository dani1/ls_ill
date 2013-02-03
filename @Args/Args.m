classdef Args < dynamicprops
% This class serves to interpret optional orguments for a function in a structured manner

methods

    function args = Args( varargin )

    cellargs = varargin;                % get the cell array

    l = length(cellargs);

    try     assert( l ~= 0 & ~mod( l, 2 ) );
    catch err;  error('Optional args not even!');
    end

    for i = 1 : floor(l/2)
        name  = cellargs{2*i -1};
        value = cellargs{2*i};

        args.addprop(name);
        args.(name) = value;
    end

    end

end

end
