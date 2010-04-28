%==========================================================================
% here the dilution is calculated starting from the pipette types and
% quantities. The input syntax is the following.
%
% 1. If the sample is diluted ONCE, just input two arrays for the types and
% values, plus the number of sample pipettes, es:
%
% dilution ( {'blue' 'yellow' 'yellow'}, [ 500 100 100 ], 1 )
%
% In this case, the first pipette is considered to be 500mul of
% blue-pipetted protein solution, the other two are water pipettes.
%
% 2. If the sample is diluted various times, just repeat the syntax ( this
% is understood through 'varargin', e.g.:
%
% dilution ( {'b' 'y'}, [ 300 100 ], 1, {'y' 'y'} [ 100 100 ], 1 )
%==========================================================================
function [ gamma dgamma ] = dilution ( varargin )

 gamma = 1; dgamma = 0;
 l = length(varargin);

 if mod(l,3)
  error('Please respect input syntax: values,names,num,... ');
 else
  for j = 1 : l/3	

   p = []; pn = []; dp = [];
   w = []; wn = []; dw = [];
   num	= varargin{3*j};
   % fill the protein
   for i = 1 : num
    pn{i}	= varargin{3*j-1}{i};
    p(i)	= varargin{3*j-2}(i);
    dp(i)	= pipetting_error(p(i),pn{i});
   end;
   % fill the water
   for i = num+1 : length(varargin{3*j-2})
    wn{i-num}	= varargin{3*j-1}{i};
    w(i-num)	= varargin{3*j-2}(i);
    dw(i-num)	= pipetting_error(w(i-num),wn{i-num});
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
