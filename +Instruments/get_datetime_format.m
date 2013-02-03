function [ regexpstr ] = get_datetime_format( str )
    % custom definition of used date formats -> for now expandable via try ... catch for consistency reasons -> in future eventually via another faster method
    % NB: use first the most common format for faster computation.
    % TODO: expand for user definitions
    try
        regexpstr =  '"mm/dd/yyyy" "HH:MM:SS PM"';
        datenum(str, regexpstr);
    catch err
        try 
            regexpstr =  '"dd/mm/yy" "HH:MM:SS"';
            datenum(str, regexpstr);
        catch err1
            try 
                regexpstr =  '"dd.mm.yy" "HH:MM:SS"';
                datenum(str, regexpstr);
            catch err2
                regexpstr = '';
        end
    end
end
