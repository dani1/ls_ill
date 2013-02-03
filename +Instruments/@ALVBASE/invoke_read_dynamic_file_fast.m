function point = invoke_read_dynamic_file_fast(self, path )
    % launch read_dynamic_file_fast( written in c), and save results in DLS.Point class.
    %--------------------------------------------------------------------------
    % change home directory to full path, since fopen does not recognize
    %it in C
    %--------------------------------------------------------------------------
    if path(1) == '~'
        if ispc %not tested
            homepath = winqueryreg('HKEY_CURRENT_USER',...
                ['Software\Microsoft\Windows\CurrentVersion\' ...
                'Explorer\Shell Folders'],'Personal');
        else
            homepath = getenv('HOME');
            path     = [homepath path(2:end)];
        end
    end
    %==========================================================================
    % get data from dynamic file
    %==========================================================================
    
    [tau g dg angle T datetime] = self.read_dynamic_file_fast( path );
    
    %==========================================================================
    % save data in DLS.Point class and correct correlation function
    %==========================================================================
    point              = DLS.Point;
    point.Instrument   = self;
    point.T            = T;
    point.Angle        = angle;
    point.Tau_raw      = tau;
    point.G_raw        = g;
    point.dG_raw       = dg;
    point.datetime_raw = datetime;
    point.correct_G();
end
