%******************************************************************************
% Written By: Daniel Soraruf
% last change : 28/04/2012
%******************************************************************************
%==============================================================================
% raw SLS data at one angle, calculate KcR
%==============================================================================
classdef AngleData < handle
	properties
		scatt_angle;
		% allocate count to 0x0 array of struct
		count = struct([]);
		mean_count_rate;
        error_mean_count_rate;
		mean_monitor_intensity;
		error_mean_monitor_intensity;
        mean_temperature;
		KcR;
		dKcR;
	end
	methods
		%----------------------------------------------------------------------
		% constructor : set angle
		%----------------------------------------------------------------------
		function self = AngleData(scatt_angle)
			self.scatt_angle = scatt_angle;
		end
		%----------------------------------------------------------------------
		% add struct of one count to class struct (members scatt_angle,
		% count_rate, minitor_intensity, error_count_rate);
		%----------------------------------------------------------------------
		function add(self, count_struct)
			len = length(self.count);
			if count_struct.scatt_angle == self.scatt_angle;
				self.count(len + 1).count_rate = count_struct.count_rate;
				self.count(len + 1).monitor_intensity = count_struct.monitor_intensity;
				self.count(len + 1).error_count_rate = count_struct.error_count_rate;
				self.count(len + 1).file_index = count_struct.file_index;
				self.count(len + 1).temperature = count_struct.temperature;
			end
		end

		% calculate mean count rate and monitor intensity with errors (std_dev)
		function e = calc_error(self)
			e = 0;
			for i = 1 : length(self.count)
				e = e + self.count(i).error_count_rate ^ 2;
			end
			e = sqrt(e) / length(self.count);
			
		end
		function calc_mean(self)
			%self.check_cr();
			if ~isempty(self.count)
				len = length(self.count);
				% set values to zero
				mean_cr = 0;
				mean_i_mon = 0;
				error_mean_cr = 0;
				error_mean_i_mon = 0;
				% calculate mean values
				for i = 1 : len
					mean_cr = mean_cr + self.count(i).count_rate;
					mean_i_mon = mean_i_mon + self.count(i).monitor_intensity;
				end
				mean_cr = mean_cr / len;
				mean_i_mon = mean_i_mon / len;
				%--------------------------------------------------------------
				% calculate errors (standard deviation)
				% TODO: check errors, since they do not match ALV ones!!!
				%--------------------------------------------------------------
				for i = 1 : len
					dev_cr = self.count(i).count_rate - mean_cr;
					error_mean_cr = error_mean_cr + dev_cr * dev_cr;
					dev_i_mon = self.count(i).monitor_intensity - mean_i_mon;
					error_mean_i_mon = error_mean_i_mon + dev_i_mon * dev_i_mon;
					
				end
				self.mean_monitor_intensity = mean_i_mon;
				self.mean_count_rate = mean_cr;
				self.error_mean_count_rate = sqrt(error_mean_cr / (len -1));
				self.error_mean_monitor_intensity = sqrt(error_mean_i_mon / (len - 1));
				self.mean_temperature = mean([self.count.temperature]);
			else
				disp(['no data at angle' self.scatt_angle])
				self.mean_count_rate = 0;
				self.mean_monitor_intensity = 0;
			end
		end
		%----------------------------------------------------------------------
		% Error propagation function for dR = 1/R !!!
		%----------------------------------------------------------------------
		function dR = R_error_propagation(self, standard, solvent, R_solution, dR_solution, angle_index)
			R_tol = standard.ratio(angle_index);
			dR_tol = standard.error_ratio(angle_index)* R_tol / 100;
			R_solv = solvent.ratio(angle_index);
			dR_solv = solvent.error_ratio(angle_index)* R_solv / 100;
			RR = standard.rayleigh_ratio(angle_index);

			% calculate error of 1 over R using gaussian error propagation
			dR =sqrt((dR_tol / (R_solution - R_solv))^2 + (R_tol * dR_solv /...
			   	(R_solution -R_solv)^2)^2 + (R_tol * dR_solution /...
			   	(R_solution - R_solv)^2)^2)/RR;
		end
		%======================================================================
		% calculate Kc/R at certain angle
		% protein_conc in mg/ml !!!
		% dn_over_dc normally in [0.1:0.3] ml/g
		% TODO: error propagation
		%======================================================================
		function calc_kc_over_r(self, standard, solvent, protein_conc, dn_over_dc, instrument)
			self.calc_mean();
            angle_tolerance = 1e-5;
            angle_index = find(abs(standard.scatt_angle-self.scatt_angle) < angle_tolerance);
			protein_conc = protein_conc * 1e-3;
			wavelength = instrument.Lambda * 1e-8;% A to cm
			%wavelength = 0.00006328 ; % cm
			number_avogadro = Constants.Na; 
			
			% calculate optical constant : using correction for cylindrical cuvettes
			% multiplication by n_std ^2 / n_solv ^2
			K = (2 * pi * dn_over_dc * standard.refraction_index )^2 /...
				(wavelength^4 * number_avogadro);
			%sample_excess_count_rate = (count_rate - solvent.count_rate(angle_index)) ...
			%/ standard.count_rate(angle_index);

			% calculate ratio of solution
			R_solution =  self.mean_count_rate * sin( self.scatt_angle * pi/180)...
			   	/ self.mean_monitor_intensity;
			dR_solution = sqrt((self.error_mean_count_rate / self.mean_monitor_intensity...
				)^2 + (self.mean_count_rate * self.error_mean_monitor_intensity /...
				self.mean_monitor_intensity^2)^2)*sin(self.scatt_angle * pi/180);

			% R needed for Kc/R, calculated using ratios saved in solvent
			% and standard file
			R = (R_solution - solvent.ratio(angle_index) ) / standard.ratio(angle_index)...
			   	* standard.rayleigh_ratio(angle_index);
			dR = self.R_error_propagation(standard, solvent, R_solution, dR_solution, angle_index);
			one_over_R = 1 / R;
			rel_error = dR / one_over_R;

			% calculate KcR and relative error
			% TODO: calculate errors and check standard deviations
			self.KcR = K * protein_conc/ R;
			self.dKcR = K * protein_conc * dR;
			%disp(self.dKcR / self.KcR)
			
		end
		%----------------------------------------------------------------------
		% show angle, mean count rate and single count rates for checks
		% TODO: rewrite for new save method. not necessary
		%----------------------------------------------------------------------
		function show_cr(self)
			a = self.scatt_angle;
			cr = [self.count.count_rate];
			mean_cr = self.mean_count_rate;
			str =[ num2str(a,'% 3.i') ' : ' num2str(mean_cr, '%.3f') ' -> '];
			for i = 1 : length(cr)
				str = [ str num2str(cr(i),'%.2f') ' '];
			end
			disp(str);
		end
		%----------------------------------------------------------------------
		% check that count_rate does not have strange scattering (e.g.dust)
		% TODO: rewrite for new save technique. !!! not explicitly needed
		%----------------------------------------------------------------------
		function check_cr(self)
			cr_tolerance_limit = 0.05;
			cr = self.count_rates;
			mean_cr = mean(cr);
			for i = 1 : length(cr)
				if abs(cr(i) - mean_cr)/mean_cr > cr_tolerance_limit
					cr_tmp = cr(:);
					cr_tmp(i) = [];
					if mean(cr_tmp) < mean_cr
						disp(['Ignored Data at ' num2str((self.scatt_angle)) ...
							': CR=' num2str(cr(i),'%.2f') ' mean=' num2str(mean_cr,'%.2f')...
							' mean_new=' num2str(cr_tmp, '%.2f')])
						self.ignore = [self.ignore ; i];
					end
				end
			end
		end
	end
end
