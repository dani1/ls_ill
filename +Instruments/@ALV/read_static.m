function [sls_point RawData] = read_static(path_standard, path_solvent, path_file, protein_conc, dn_over_dc, start_index, end_index, count_number, varargin)
	% input: (path_standard, path_solvent, path_file, protein_conc, dn_over_dc, start_index {autosave file} , end_index {autosave file})
	% 
	% read and calculate data for static light scattering, giving as input
	% standard and water .tol files, as well as the autosave path (without dddd.ASC)
	% and protein concentration, differential index of refraction. 
	% furthermore start_index and end_index are the file indices which to load
	% 
	% output : [Point SlsData]
	% SlsData is AngleData class array with single file-data info

% ******************************************************************************
% Written By: Daniel Soraruf
% last change : 19/02/2012
% ******************************************************************************
	% for debug purpose commented variable definitions
	% path_solvent = '~/Documents/tesi/data/data_raw/LS/2011_11_04/Water.tol';
	% path_standard = '~/Documents/tesi/data/data_raw/LS/2011_11_04/Toluene.tol';
	% path_file = '~/Documents/tesi/data/data_raw/LS/2011_11_04/BSA_1gl_NaCl_200mM';
	% protein_conc = 0.001;
	% dn_over_dc = 0.1;
	% start_index = 0;
	% end_index = 37;
	% ------------------------------------------------------------------------------
	% get solvent and standard data
	% ------------------------------------------------------------------------------
	solvent = Instruments.ALV.read_tol_file(path_solvent);
	standard = Instruments.ALV.read_tol_file(path_standard);
	index = 0;
	point(end_index - start_index + 1) = struct('scatt_angle', 0, 'count_rate', 0, 'monitor_intensity', 0, 'error_count_rate', 0);
	% ------------------------------------------------------------------------------
	% get data from autosave ALV files
	% ------------------------------------------------------------------------------
	if path_file(1) == '~'
		homepath = getenv('HOME'); % !!! works only on unix systems
		path_file = [homepath path_file(2:end)];
		% disp(path);
		% path = '/Users/daniel/Documents/tesi/data/data_raw/LS/2011_11_04/BSA_1gl_NaCl_200mM0003.ASC';
	end
	for i = start_index : end_index
		index = index + 1;
		file = [ path_file num2str(i,'%4.4u') '.ASC' ];
		%[count_rate1 count_rate2 point(index).monitor_intensity point(index).scatt_angle point(index).temperature datetime]...
		% = Instruments.ALV.read_static_from_autosave_fast(file);
        [count_rate1 count_rate2 point(index).monitor_intensity point(index).scatt_angle point(index).temperature datetime]...
		= Instruments.ALV.read_static_from_autosave(file);
		point(index).cr1 = count_rate1;
		point(index).cr2 = count_rate2;
		point(index).count_rate = count_rate1 + count_rate2;
		point(index).error_count_rate = sqrt(count_rate1 * 1000) + sqrt(count_rate2 * 1000);
		point(index).file_index = index;
		point(index).datetime = datenum(datetime);
	end

%	[KcR] = calc_kc_over_r(scatt_angle, standard, solvent,cr_mean,0.001, Imean);

	% find all angles in file
	a = unique([point.scatt_angle]);
%	sort them ( to be consistent with the program generated file)
	angles = sort(a);
	% preallocate memory
	SlsData = SLS.AngleData.empty(length(angles),0);
	angle_tolerance = 1e-3;
	for index = 1 : length(angles)
		
		SlsData(index) = SLS.AngleData(angles(index));
		for index_1 = 1 : length(point);
			% save data at one ANGLE in SlsData
			if abs(point(index_1).scatt_angle - angles(index)) < angle_tolerance
				SlsData(index).add(point(index_1));
			end
		end
	end
	% --------------------------------------------------------------------------
	% calc Kc over R
	% --------------------------------------------------------------------------
	for i = 1 : length(SlsData)
		SlsData(i).calc_kc_over_r(standard, solvent, protein_conc, dn_over_dc,Instruments.ALV);
		% SlsData(i).show_cr();
	end

	% --------------------------------------------------------------------------
	% Save into SLS.Point array
	% --------------------------------------------------------------------------
	% disp(SlsData(1).KcR)
	for i = 1 : length(SlsData)
		sls_point(i) = SLS.Point;
		sls_point(i).Instrument = Instruments.ALV;
		sls_point(i).T = SlsData(i).mean_temperature;
		sls_point(i).Angle = SlsData(i).scatt_angle;
		sls_point(i).KcR_raw = SlsData(i).KcR;
		sls_point(i).dKcR_raw = SlsData(i).dKcR;
		sls_point(i).datetime = mean([SlsData(i).count(1:end).datetime]);
    end
	RawData.SlsData = SlsData;
	RawData.solvent = solvent;
	RawData.standard = standard;
end
