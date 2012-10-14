function [ s e nc] = find_start_end ( path )
%out:[s e nc] -> s = start, e = end , nc = number of counts for each angle
    i = -1;
    flag = 0;
    while flag == 0 && i < 200
        i = i + 1;
        flag = exist( [ path num2str(i,'%4.4u') '_0001' '.ASC'] );
    end
    if flag == 0
        error(['DLS files: "' path '" not found']);
        s = -1;
        return
    else
        s = i;
    end
    flag = 1;
    while flag ~= 0 && i < 500
        i = i + 1;
        flag = exist( [ path num2str(i, '%4.4u') '_0001' '.ASC' ]);
    end
    e = i - 1;
    flag = 1;
    i = 0;
    while flag ~= 0 && i < 3000
        i = i + 1;
        flag = exist( [ path num2str(s, '%4.4u') '_' num2str(i, '%4.4u') '.ASC' ]);
    end
    nc = i - 1;
    % check whether every angle has same amount of counts 
    index_of_different_count_numbers = 1;
    nc_array(1) = nc;
    s_array(1) = s;
    for i_a = s : e
        i = 0;
        flag = 1;
        while flag ~= 0 
            i = i + 1;
            flag = exist( [ path num2str(i_a, '%4.4u') '_' num2str(i, '%4.4u') '.ASC' ]);
        end
        if i~=nc_array(index_of_different_count_numbers)+1
            if (i_a < e)
                s_array(index_of_different_count_numbers + 1) = i_a;
            end
            e_array(index_of_different_count_numbers) = i_a - 1;
            nc_array(index_of_different_count_numbers + 1) = i - 1;
            index_of_different_count_numbers = index_of_different_count_numbers + 1;
            % exit for cycle
        end
        e_array(index_of_different_count_numbers) = i_a;
    end
    s = s_array;
    e = e_array;
    nc = nc_array;
end
