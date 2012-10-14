function s = read_tol_file(path_of_tol_file)
%        path_of_tol_file = '~/Documents/tesi/data/data_raw/LS/2011_10_31/Water2.tol'
    %path_of_tol_file
    tol = importdata(path_of_tol_file,'\t', 3);
    s.scatt_angle          = tol.data(:,1); % in degrees
    s.q2_scatt             = tol.data(:,2); % in 1/m^2
    s.count_rate           = tol.data(:,3); % in percent
    s.error_count_rate     = tol.data(:,4); % in percent
    s.temperature          = tol.data(:,7); % in K
    s.error_temperature    = tol.data(:,8); % in percent
    s.ratio                = tol.data(:,9);
    s.error_ratio          = tol.data(:,10);
    s.rayleigh_ratio       = tol.data(:,11); % in 1/cm
    s.error_rayleigh_ratio = tol.data(:,12); % in percent
    s.refraction_index     = str2num(tol.textdata{3});
end
