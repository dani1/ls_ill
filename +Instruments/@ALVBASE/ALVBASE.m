classdef ALVBASE
% This class describes the features of the ALV light scattering instrument at ILL

 properties ( Constant )

 end
 methods

  Point = invoke_read_dynamic_file_fast(self, path);
  Point	= read_dynamic_file	(self, path );
  Point	= read_static_file 	(self, path );

 end

 methods ( Static )

  % [Point RawData] = read_static(path_standard, path_solvent, path_file, protein_conc, dn_over_dc, start_index, end_index, count_number);
  [t gt dgt Angle temperature datetime] = read_dynamic_file_fast( path );
  s = read_tol_file(path_of_tol_file);
  [count_rate1 count_rate2 I_mon angle temperature datetime] = read_static_from_autosave(path_of_autosave_file);
  [count_rate1 count_rate2 I_mon angle temperature datetime] = read_static_from_autosave_fast(path_of_autosave_file);
  % [s e nc] = find_start_end( path );
  % [fname] = generate_filename( path_file,  angle_index, count_index);
 end

end
