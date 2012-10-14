function fname = generate_filename(path_file, index, varargin)
    fname = [path_file num2str(index, '%4.4u') '.ASC'];
