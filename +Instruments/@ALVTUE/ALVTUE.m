classdef ALVTUE < Instruments.ALVBASE
% This class describes the features of the ALV light scattering instrument at IAP Tübingen
properties
	attenuator;
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

%  methods
% function self = ALVTUE ( self )
% 	self.attenuator = self.get_attenuator_corrections();
%  end
%  end
 methods ( Static )

	 % read static from table
  % [Point RawData] = read_static(path_standard, path_solvent, path_file, protein_conc, dn_over_dc, start_index, end_index, count_number, varargin)
  % read static files from autosave
  % s = read_tol_file(path_of_tol_file)
  % [count_rate1 count_rate2 I_mon angle temperature datetime] = read_static_from_autosave(path_of_autosave_file)
  % [count_rate1 count_rate2 I_mon angle temperature datetime] = read_static_from_autosave_fast(path_of_autosave_file)
  % read dynamic files
  % Point	= read_dynamic_file	( path )
  % Point	= invoke_read_dynamic_file_fast( path )
  % [t gt dgt Angle temperature datetime] = read_dynamic_file_fast( path );
  [ s e nc] = find_start_end ( path )
  fname = generate_filename(path_file, angle_index, count_index)
function [att] = get_attenuator_corrections()
	att(7).monitor_intensity = 2960000;
	att(7).intensity_correction = 1;
	att(7).percent_transmission = 100;
	
	att(6).monitor_intensity = 870000;
	att(6).intensity_correction = 1.06   * att(7).intensity_correction;
	att(6).percent_transmission = 33;

	att(5).monitor_intensity = 315000;
	att(5).intensity_correction = 1.039   * att(6).intensity_correction;
	att(5).percent_transmission = 10;

	att(4).monitor_intensity = 120657;
	att(4).intensity_correction = 1.0713 * att(5).intensity_correction;
	att(4).percent_transmission = 3;

	att(3).monitor_intensity = 31200;
	att(3).intensity_correction = 1.169    * att(4).intensity_correction;
	att(3).percent_transmission = 1;

	att(2).monitor_intensity = 5700;
	att(2).intensity_correction = 1.8      * att(3).intensity_correction;
	att(2).percent_transmission = 0.3;

	att(1).monitor_intensity = 1;
	att(1).intensity_correction = 1.06   * att(2).intensity_correction;
	att(1).percent_transmission = 0.1;

end
 end

end
