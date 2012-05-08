classdef ALVTUE
% This class describes the features of the ALV light scattering instrument at ILL

properties 
  attenuator ;
end
 properties ( Constant )

  Goniometer	= 'ALV-CGS3';
  Correlator	= 'ALV-7004/FAST';
  Lambda	= 6328;
  Unit_Lambda	= 'A';
  Angles	= struct(	'Min',	12,	...
				'Max',	150	);
  T		= struct(	'Min',	-Constants.T0,	...	% Minimal T = 0°C
				'Max',	-Constants.T0+45	);	% Maximal T = 45°C
%attenuators values	
 end

 methods
function self = ALVTUE ( self )
	self.attenuator(1).monitor_intensity = 1;
	self.attenuator(1).intensity_correction = 1.06;
	self.attenuator(1).percent_transmission = 0.1;

	self.attenuator(2).monitor_intensity = 5700;
	self.attenuator(2).intensity_correction = 1.06;
	self.attenuator(2).percent_transmission = 0.3;

	self.attenuator(3).monitor_intensity = 31200;
	self.attenuator(3).intensity_correction = 1.06;
	self.attenuator(3).percent_transmission = 1;

	self.attenuator(4).monitor_intensity = 120657;
	self.attenuator(4).intensity_correction = 1.15;
	self.attenuator(4).percent_transmission = 3;

	self.attenuator(5).monitor_intensity = 315000;
	self.attenuator(5).intensity_correction = 1.1;
	self.attenuator(5).percent_transmission = 10;

	self.attenuator(6).monitor_intensity = 870000;
	self.attenuator(6).intensity_correction = 1.06;
	self.attenuator(6).percent_transmission = 33;

	self.attenuator(7).monitor_intensity = 2960000;
	self.attenuator(7).intensity_correction = 1;
	self.attenuator(7).percent_transmission = 100;
 end
 end
 methods ( Static )

	 % read static from table
  Point	= read_static_file 	( path )
  % read static files from autosave
  [Point SlsData RawData] = read_static(path_standard, path_solvent, path_file, protein_conc, dn_over_dc, start_index, end_index, count_number)
  s = read_tol_file(path_of_tol_file)
  [count_rate1 count_rate2 I_mon angle temperature] = read_static_from_autosave(path_of_autosave_file)
  [count_rate1 count_rate2 I_mon angle temperature datetime] = read_static_from_autosave_fast(path_of_autosave_file)
  % read dynamic files
  Point	= read_dynamic_file	( path )
  Point	= invoke_read_dynamic_file_fast( path )
  [t gt dgt Angle temperature] = read_dynamic_file_fast( path );
  [ s e nc] = find_start_end ( path )
  
 end

end
