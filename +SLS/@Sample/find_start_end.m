function [ s e ] = find_start_end ( path )

 i = -1;
 flag = 0;
 while flag == 0 & i < 200
  i	= i+1;
  flag	= exist( [ path num2str(i,'%4.4u') '.ASC' ] );
 end

 if flag == 0
  error('DLS files not found');
 else
  s	= i;
 end

 flag = 1;
 while flag ~= 0 & i < 500
  i	= i+1;
  flag	= exist( [ path num2str(i,'%4.4u') '.ASC' ] );
 end

 e	= i-1;

end
