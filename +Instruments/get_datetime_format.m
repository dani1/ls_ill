function [ regexpstr ] = get_datetime_format( str )
    try
        regexpstr =  '"dd/mm/yyyy" "HH:MM:SS PM"';
        datenum(str, regexpstr);
    catch err
        try 
            regexpstr =  '"dd.mm.yy" "HH:MM:SS PM"';
            datenum(str, regexpstr);
        catch err1
            try 
                regexpstr =  '"dd-mm-yy" "HH:MM:SS PM"';
                datenum(str, regexpstr);
            catch err2
                regexpstr = '';
        end
    end
end

