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
  raw_data_path
  date_experiment
  RawData

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
 if any(strcmp('path_standard', properties(a)))
	   bool_get_data_from_autosave = 1;
else
	   bool_get_data_from_autosave = 0;
	end
	% get data from table
	if ~bool_get_data_from_autosave
	   try		self.Point	= self.Instrument.read_static_file ( a.Path );  	% get the KcR and angles
	   catch disp(self)
	       error('Error loading the static file!');
	    
	   end
    else
        disp(['Load SLS:' a.Path])
		[start_index, end_index, nc] = self.Instrument.find_start_end( a.Path );
		path_standard = a.path_standard;
		path_solvent  = a.path_solvent;
		[self.Point self.RawData] = self.Instrument.read_static(path_standard, path_solvent, ...
			a.Path, self.C, self.dndc, start_index, end_index, nc);
	end
   self.raw_data_path = a.Path;
   pointprops	= {	'Protein',	'Salt',		...
			'C',		'C_set',	...
			'Cs',				...
			'n',		'n_set',	...
			'dndc',		'dndc_set'	};
   for i = 1 : length(pointprops)
     [ self.Point.(pointprops{i}) ]	= deal(a.(pointprops{i}));			% the deal function rocks!
   end

  end

  function correct_attenuator( self )
	  if any(strcmp('SlsData', properties(self.RawData)))
		sls_data = RawData.SlsData;
		for i_angle = 1 : length(sls_data)
			for i_att = 1 : length(self.Instrument.attenuator)
				if ( round(sls_data(i_angle).count(1).monitor_intensity / self.Instrument.attenuator(i_att).monitor_intensity) == 1)
					sls_data(i_angle).KcR = sls_data.KcR * self.Instrument.attenuator(i_att).intensity_correction;
					sls_data(i_angle).dKcR = sls_data.dKcR * self.Instrument.attenuator(i_att).intensity_correction;
					% exit i_att cycle
					break
				end
			end
			self.Point(i_angle).KcR = sls_data(i_angle).KcR
			self.Point(i_angle).dKcR = sls_data(i_angle).dKcR
		end
	  else
		  disp('no raw data from autosave available ! -> exit')
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
 methods(Access = private, Static)
	 [s e nc] = find_start_end( path )
 end

end
