%==========================================================================
% here the dilution is calculated starting from the pipette types and
% quantities. The input syntax is the following.
%
% 1. If the sample is diluted ONCE, just input two arrays for the types and
% values, plus the number of sample pipettes, es:
%
% dilution ( [500 100 100], {'blue' 'yellow' 'yellow'}, 1 )
%
% In this case, the first pipette is considered to be 500mul of
% blue-pipetted protein solution, the other two are water pipettes.
%
% 2. If the sample is diluted various times, just repeat the syntax ( this
% is understood through 'varargin', e.g.:
%
% dilution ( [300 100], {'b' 'y'}, 1, [100 100], {'y' 'y'}, 1 )
%
% In the future, it will be possible to input a struct instead of an array,
% but this is work in progress.
%==========================================================================
function [ gamma dgamma ] = dilution ( varargin )

 gamma = 1; dgamma = 0;

 % this is good practice because the string "varargin" should be only used as input argument of the function
 argvect = varargin;

 try									% recognition of structs as input (temporary version)
  if isstruct(argvect{:})
   fprintf('Hey, this is a struct!\n');
  end
 end

 if mod(length(argvect),3)						% if the input syntax is not correct, tell the user
  error('Please respect input syntax: values,names,num,... ');
 else

  l = length(argvect)/3;						% this is the number of dilution steps

  for j = 1 : l								% repeat the core function for every dilution step

   starting_index = 3*(j-1);						% this is the starting index for every dilut step

   p = []; pn = []; dp = [];						% empty the temporary variables
   w = []; wn = []; dw = [];						% empty the temporary variables

   num_p	= argvect{3*j};						% this is the number of sample pipettes
   num_w	= length(argvect{ starting_index+1 }) - num_p;		% this is the number of water pipettes

   % fill the protein (pipettes from 1 up to num)
   for i = 1 : num_p
    p(i)	= argvect{ starting_index+1 }(i);			% the first ones are always the pipettes
    pn{i}	= argvect{ starting_index+2 }{i};			% the second ones are the colours
    dp(i)	= pipetting_error(p(i),pn{i});				% calculate the errors
   end;

   % fill the water (remaining pipettes)
   for i = num_p+1 : num_p+num_w
    w(i-num_p)	= argvect{ starting_index+1 }(i);			% the first ones are always the pipettes
    wn{i-num_p}	= argvect{ starting_index+2 }{i};			% the second ones are the colours
    dw(i-num_p)	= pipetting_error(w(i-num_p),wn{i-num_p});		% calculate the errors
   end

   % calculate the new dilution every time, using tmp vars
   gamma_old	= gamma;
   dgamma_old	= dgamma;
   gamma_new	= sum(p) / sum( [ p w ]);
   dgamma_new	= gamma_new * sqrt(	sum( ( [ dp dw ] ./ [ p w ] ).^2 ) );
   gamma	= gamma_new * gamma_old;
   dgamma	= gamma * sqrt( ( dgamma_new / gamma_new ) ^2 + ( dgamma_old / gamma_old ) ^2 );

  end
 end

end	% dilution function
