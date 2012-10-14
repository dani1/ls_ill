function [s e nc] = find_start_end ( path )
    [s e nc] = find_start_end1 (path);
    if s < 0 
        [s e nc] = find_start_end2(path);
    end
    if s < 0
        error('DLS files not found');
    end
end
function [ s e nc] = find_start_end1 ( path )
 nc = -1;
 e = -1;
 i = -1;
 flag = 0;
 while flag == 0 && i < 200
  i	= i+1;
  flag	= exist( [ path num2str(i,'%4.4u') '.ASC' ] );
 end

 if flag == 0
  %error('DLS files not found');
  s = -1;
  return
 else
  s	= i;
 end

 flag = 1;
 while flag ~= 0 && i < 500
  i	= i+1;
  flag	= exist( [ path num2str(i,'%4.4u') '.ASC' ] );
 end

 e	= i-1;

end
function [ s e nc] = find_start_end2 ( path )
%out:[s e nc] -> s = start, e = end , nc = number of counts for each angle
    i = -1;
    flag = 0;
    while flag == 0 && i < 200
        i = i + 1;
        flag = exist( [ path num2str(i,'%4.4u') '_0001' '.ASC'] );
    end
    if flag == 0
        %error('DLS files not found');
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
    while flag ~= 0 && i < 200
        i = i + 1;
        flag = exist( [ path num2str(s, '%4.4u') '_' num2str(i, '%4.4u') '.ASC' ]);
    end
    nc = i - 1;
    i_a = s;
    flag = 1;
    % check whether every angle has same amount of counts 
    for i_a = s : e
        i = 1;
        while flag ~= 0 && i <= nc
            flag = exist( [ path num2str(i_a, '%4.4u') '_' num2str(i, '%4.4u') '.ASC' ]);
            i = i + 1;
        end
        if flag == 0
            e = i_a - 1;
            % exit for cycle
            break;
        end
    end
end
