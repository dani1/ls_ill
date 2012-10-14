classdef Sample
% This class takes as input some info about a sample and stores it.

properties

    Protein         % what protein?
    Salt            % what salt?
    C               % protein concentration
    Unit_C = 'g/l';
    Cs              % salt concentration
    Unit_Cs = 'mM'

    T
    Unit_T = 'K';
    n

    Point           % datapoints

end

properties ( SetAccess = private, Dependent )

    Angle
    Q

end

properties ( Hidden)

    raw_data_path
    Instrument
    C_set
    n_set
    datetime
    number_of_counts
    start_index
    end_index

end

methods

    %constructor
    function self = Sample( varargin )
        a = Args(varargin{:});% get  the args
        try self.Instrument = Instruments.(a.Instrument); % get the instrument
        catch err; error('Instrument not found!');
        end

        props = {	'Protein',	'Salt',		...
                'C' , 'C_set', ...
                'Cs', 'T'    , ...
                'n' ,  'n_set'  };
        for i = 1 : length(props)
            try		self.(props{i})		= a.(props{i});
            catch ;	warning(['Property ' props{i} ' not found!']);
            end
        end

        self.raw_data_path = a.Path;
        if any(strcmp('filegroup_index', properties(a)))
            filegroup_index = a.filegroup_index;
        else
            filegroup_index = 1;
        end
        [s_array e_array nc_array] = self.Instrument.find_start_end( self.raw_data_path );
        s = s_array(filegroup_index);
        e = e_array(filegroup_index);
        nc = nc_array(filegroup_index);
        if any(strcmp('start_index', properties(a))) && a.start_index > 0
            s = a.start_index;
        end
        if any(strcmp('end_index', properties(a))) &&a.end_index > 0
            e = a.end_index;
        end
        if any(strcmp('number_of_counts', properties(a))) && a.number_of_counts > 0
            nc = a.number_of_counts;
        end
        self.Point = DLS.Point;
        disp(['load: ' self.raw_data_path '[' num2str(s, '%4.4u') ':' num2str(e, '%4.4u') ']' ]);
        self.Point = DLS.Point;
        counter = 0;
        for i = s : e
            flag = true;
            i_c = 1;
            while flag
                counter = counter + 1;
                file = self.Instrument.generate_filename(self.raw_data_path, i, i_c);
                self.Point(counter)	= self.Instrument.invoke_read_dynamic_file_fast( file );
                % self.Point(i-s+1)	= self.Instrument.read_dynamic_file(file);
                % condition here used to simulate behavior of do - while loop of C.
                if i_c >= nc
                    flag = false;
                else
                    i_c = i_c + 1;
                end
            end
        end
        self.start_index = s;
        self.end_index = e;
        self.number_of_counts = nc;
        
        regexpstr = Instruments.get_datetime_format(self.Point(1).datetime_raw);
        if ~isempty(regexpstr)
            for i = 1 : length(self.Point)
                self.Point(i).datetime = datenum(self.Point(i).datetime_raw, regexpstr);
            end
        end
        pointprops	= {'Protein','Salt',...
                    'C' , 'C_set'    , ...
                    'Cs', ...
                    'n' , 'n_set'  };
        for i = 1 : length(pointprops)
            [ self.Point.(pointprops{i}) ]	= deal(self.(pointprops{i})); % the deal function rocks!
        end
        self.datetime = mean(horzcat(self.Point(1:end).datetime));
    end

    function Angle= get.Angle ( self )
        Angle= unique([self.Point.Angle]);
    end

    function Q = get.Q ( self )
        Q = unique([self.Point.Q]);
    end
    function Q = Qv ( self )
    % function which return full Q vector of the length of self.Point
        Q = [self.Point.Q];
    end
    function [fit_val, error_fit_val] = get_fit(self, method, parameter)
        % get_fit : function to retrieve fit values and errors of 95% confidence interval
        % input : method (e.g. 'DoubleBKG') , parameter (e.g. 'Gamma1')
        fitmethod = ['Fit_' method];
        varnames  = coeffnames(self.Point(1).(fitmethod));
        ind       = strcmp(varnames, parameter);
        index     = find(ind, 1);
        i_err     = 2 * index - 1;
        for i = 1 : length(self.Point)
            p = self.Point(i).(fitmethod);
            fit_val(i)       = p.(parameter);
            errors           = confint(p);
            error_fit_val(i) = abs((errors(i_err+1) - errors(i_err))) * 0.5;
            if isnan(error_fit_val(i))
                disp(['confidence value not defined at: ' num2str(i)])
                error_fit_val(i) = fit_val(i);
            end
        end
    end
end

methods ( Access = private, Static )

    [ s e nc] = find_start_end ( path );

end

 % FIT METHODS
methods

    function fit ( self , model )
        for i = 1 : length( self.Point )
            self.Point(i).fit( model );
        end
    end
    function fit_raw ( self , model )
        for i = 1 : length( self.Point )
            self.Point(i).fit_raw( model );
        end
    end
    function invert_laplace ( self )
        for i = 1 : length( self.Point )
            fprintf([num2str(i) ': ']);
            self.Point(i).invert_laplace;
            fprintf('\n');
        end
    end
end
end
