classdef ExperimentData
    
    properties
        ExperimentName = '';
        Path = '';
        Protocol = [];
        % is it ok to load void strings?
        Standard = [];
        Solvent = [];
        Sample = [];
    end
    
    methods
        % constructor for class. Just load everything for now
        function ed = ExperimentData( experimentName, path, protocol, standard, solvent)
            % export parameters to the class
            ed.ExperimentName = experimentName;
            ed.Path = path;
            ed.Protocol = protocol;
            ed.Standard = standard;
            ed.Solvent = solvent;
            
            % preallocate data
            