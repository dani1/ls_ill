function nu = nuance ( color, l )
% output l nuances of the chosen color

 c		= sqrt( ( [ 1 : l ] - 0.1 * ones(1,l) ) / l );
 n0(1:l)	= 0;
 n1(1:l)	= 0.1;
 n2(1:l)	= 0.2;
 n3(1:l)	= 0.3;

 switch color
  case 'r'
   nu = [ c' n0' n0' ];

  case 'g'
   nu = [ n1' c' n2' ];

  case 'b'
   nu = [ n1' n2' c' ];

  case 'gray'
   nu = [ c' c' c' ];

 end

end
