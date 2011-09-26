function outclass = filterout(e,index_conc,indices_abscissa)
 % e is the class
 % index_conc is the figure index (concentration)
 % indices_abscissa are the point indices
 %
 % The idea is the following: I get the points from the indices, then I
 % delete the points by looking exactly for them
 outclass = e;

 indices_real = find([outclass.Point.C] == outclass.C(index_conc));	% restrict the indices to the right concentration
 indices_real = indices_real( indices_abscissa );			% restrict the indices to the right abscissa number
 point = outclass.Point( indices_real );				% get the points

 outclass.Point = outclass.Point( ~ismember(outclass.Point,point) );	% select only the points that are not member of p

end
