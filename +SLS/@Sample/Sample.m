classdef Sample
% This class takes as input some info about a sample and stores it.

 properties

  Protein					% what protein?
  Salt						% what salt?
  C						% protein concentration
  Unit_C	= 'g/l';
  Cs						% salt concentration
  Unit_Cs	= 'mM'

  Unit_KcR	= 'mol g^{-1}'
  Unit_X_T	= 'l * J^{-1}';

  T
  Unit_T	= 'K';
  n
  dndc

  Point						% datapoints

 end

 properties ( SetAccess = private, Dependent )

  Angle
  Q
  KcR
  dKcR
  X_T
  dX_T

 end

 properties ( Hidden)

  Instrument
  C_set
  n_set
  dndc_set

 end

 methods

  function self = Sample( varargin )

   a	= Args(varargin{:})	;							% get the args

   try		self.Instrument	= Instruments.(a.Instrument);				% get the instrument
   catch 	error('Instrument not found!');
   end

   props	= {	'Protein',	'Salt',		...
			'C',		'C_set',	...
			'Cs',		'T',		...
			'n',		'n_set',	...
			'dndc',		'dndc_set'	};
   for i = 1 : length(props)
    try		self.(props{i})		= a.(props{i});
    catch	warning(['Property ' props{i} ' not found!']);
    end
   end

   try		self.Point	= self.Instrument.read_static_file ( a.Path );  	% get the KcR and angles
   catch disp(self)
       error('Error loading the static file!');
      
   end

   pointprops	= {	'Protein',	'Salt',		...
			'C',		'C_set',	...
			'Cs',				...
			'n',		'n_set',	...
			'dndc',		'dndc_set'	};
   for i = 1 : length(pointprops)
     [ self.Point.(pointprops{i}) ]	= deal(a.(pointprops{i}));			% the deal function rocks!
   end

  end

  function KcR	= get.KcR ( self )
   y	= [ self.Point.KcR ];
   w	= 1./ [ self.Point.dKcR ].^2;
   KcR	= sum( y .* w ) / sum( w );
  end 

  function dKcR	= get.dKcR ( self )
   y	= [ self.Point.dKcR ];
   w	= 1./ [ self.Point.dKcR ].^2;
   dKcR	= sum( y .* w ) / sum( w );
  end

  function X_T	= get.X_T ( self )
   y	= [ self.Point.X_T ];
   w	= 1./ [ self.Point.dX_T ].^2;
   X_T	= sum( y .* w ) / sum( w );   
  end

  function dX_T	= get.dX_T ( self )
   y	= [ self.Point.dX_T ];
   w	= 1./ [ self.Point.dX_T ].^2;
   dX_T	= sum( y .* w ) / sum( w );   
  end

  function Angle= get.Angle ( self )
   Angle= unique([self.Point.Angle]);
  end

  function Q	= get.Q ( self )
   Q	= unique([self.Point.Q]);
  end

 end

end
