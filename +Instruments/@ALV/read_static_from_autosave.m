function [count_rate1 count_rate2 I_mon angle temperature] = read_static_from_autosave(path_of_autosave_file)
	% reads the mean count rate 0,1 , Monitor diode intensity and angle from autosave ASCII file
	%path_of_autosave_file =
	%'~/Documents/tesi/data/data_raw/LS/2011_10_31/BSA_5gl_0004.ASC'
	fid = fopen(path_of_autosave_file);
	str = fgetl(fid);
    while ~strcmp(str(1:11), 'Temperature')
		str = fgetl(fid);
	end
	[ tmp tmp tmp temperature ] = strread(str, '%s %s %s %f'); % read angle from line
	while ~strcmp(str(1:5), 'Angle')
		str = fgetl(fid);
	end
	[ tmp tmp tmp angle ] = strread(str, '%s %s %s %f'); % read angle from line
	while ~strcmp(str(1:7), 'MeanCR0')
		str = fgetl(fid);
	end
	[tmp tmp tmp count_rate1] = strread(str,'%s %s %s %f'); % read CR1 from line
	str = fgetl(fid);
	[tmp tmp tmp count_rate2] = strread(str,'%s %s %s %f'); % read CR2 from line
	str = fgetl(fid);
	while ~feof(fid);
		if strfind(str, 'Monitor Diode')
			[tmp tmp I_mon] = strread(str, '%s %s %f');
			break
		end
		str = fgetl(fid);
	end
	fclose(fid);
end
