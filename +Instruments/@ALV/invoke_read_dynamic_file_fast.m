function point = invoke_read_dynamic_file_fast( path )
	% launch read_dynamic_file_fast( written in c), and save results in DLS.Point class.
	%--------------------------------------------------------------------------
	% change home directory to full path, since fopen does not recognize
	%it in C
	%--------------------------------------------------------------------------
	if path(1) == '~'
		homepath = getenv('HOME'); % !!! works only on unix systems
		path = [homepath path(2:end)];
		%disp(path);
		%path = '/Users/daniel/Documents/tesi/data/data_raw/LS/2011_11_04/BSA_1gl_NaCl_200mM0003.ASC';
	end
	%==========================================================================
	% get data from dynamic file
	%==========================================================================
	
	[tau g dg angle T] = Instruments.ALV.read_dynamic_file_fast( path );
	
	%==========================================================================
	% save data in DLS.Point class and correct correlation function
	%==========================================================================
	point = DLS.Point;
	point.Instrument = Instruments.ALV;
	point.T = T;
	point.Angle = angle;
	point.Tau_raw = tau;
	point.G_raw = g;
	point.dG_raw = dg;		% this triggers the event in Point
	point.correct_G; % is function !
end
