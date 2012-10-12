function fname = generate_filename(path_file, angle_index, count_index, varargin)
	fname = [path_file num2str(angle_index, '%4.4u') '_' num2str(count_index, '%4.4u') '.ASC'];
