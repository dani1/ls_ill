classdef Sample
% This class takes as input some info about a sample and stores it.

 properties

  Protein					% what protein?
  Salt						% what salt?
  C						% protein concentration
  Unit_C	= 'g/l';
  Cs						% salt concentration
  Unit_Cs	= 'mM'

  T
  Unit_T	= 'K';
  n

  Point						% datapoints

 end

 properties ( SetAccess = private, Dependent )

  Angle
  Q

 end

 properties ( Hidden)

  Instrument
  C_set
  n_set
  raw_data_path
  date_experiment

 end

 methods

  function self = Sample( varargin )

   a	= Args(varargin{:});% get the args
   try		self.Instrument	= Instruments.(a.Instrument);				% get the instrument
   catch ; error('Instrument not found!');
   end

   props	= {	'Protein',	'Salt',		...
			'C',		'C_set',	...
			'Cs',		'T',		...
			'n',		'n_set'		};
   for i = 1 : length(props)
    try		self.(props{i})		= a.(props{i});
    catch ;	warning(['Property ' props{i} ' not found!']);
    end
   end

   self.raw_data_path = a.Path;
   [ s e nc] = self.find_start_end( a.Path );
	disp(['load: ' a.Path '[' num2str(s, '%4.4u') ':' num2str(e, '%4.4u') ']' ]);
   
   self.Point	= DLS.Point;
   if nc < 0
	   for i = s : e
		file		= [ a.Path num2str(i,'%4.4u') '.ASC' ];
		self.Point(i-s+1)	= self.Instrument.invoke_read_dynamic_file_fast( file );
		%self.Point(i-s+1)	= self.Instrument.read_dynamic_file(file);
	   end
   else
       counter = 0;
	   for i = s : e
		   for i_c = 1 : nc
            counter = counter + 1;
			file		= [ a.Path num2str(i,'%4.4u') '_' num2str(i_c, '%4.4u') '.ASC' ];
			self.Point(counter)	= self.Instrument.invoke_read_dynamic_file_fast( file );
			%self.Point(i-s+1)	= self.Instrument.read_dynamic_file(file);
			end
	   end
   end

   pointprops	= {	'Protein',	'Salt',		...
			'C',		'C_set',	...
			'Cs',				...
			'n',		'n_set'		};
   for i = 1 : length(pointprops)
     [ self.Point.(pointprops{i}) ]	= deal(a.(pointprops{i}));			% the deal function rocks!
   end

  end

  function Angle= get.Angle ( self )
   Angle= unique([self.Point.Angle]);
  end

  function Q	= get.Q ( self )
   Q	= unique([self.Point.Q]);
  end

 end

 methods ( Access = private, Static )

  [ s e nc] = find_start_end ( path );

 end

 % FIT METHODS
 methods

  function fit ( self , model )
   for i = 1 : length( self.Point )
    self.Point(i).fit( model );
   end
  end

  function invert_laplace ( self )
   for i = 1 : length( self.Point )
    fprintf([num2str(i) ': ']);
    self.Point(i).invert_laplace;
    fprintf('\n');
   end
  end

 end

end
