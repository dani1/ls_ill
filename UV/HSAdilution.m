%==========================================================================
% here the dilution is calculated starting from the pipette types and
% quantities. The input syntax is the following:
%
% For each sample one struct like the following example:

% struct('PipetteForProtein', ['y' 'y'],'volumeProtein',[120 120],...
% 'PipetteForWater',['b'], 'volumeWater', [760])

%if you have diluted several times, use one row in the field matrix for
%each dilution step:
% struct('PipetteForProtein', ['y' 'y'],'volumeProtein',[70; 100],...
% 'PipetteForWater',['b' 'b'], 'volumeWater', [350; 300])

 

function [ gamma dgamma ] = HSAdilution ( inputstruct )

 gamma = 1; dgamma = 0;
 %arrays for the dilution error
 %Protein
    p  = inputstruct.volumeProtein;         % p is a matrix containing the volume of protein with one row for each dilution step 
    pn = inputstruct.PipetteForProtein;     % pn is an array containing the pipetts used
    dp = zeros(size(p));                    % dp is the Matrix for the pipetting error
    
    for i = 1 : numel(pn)
        if (p(i)>=2)
            dp(i)	= pipetting_error(p(i),pn(i));  %the pipetting error is calculated 
        end
    end;

  %Water
    w  = inputstruct.volumeWater;
    wn = inputstruct.PipetteForWater;
    dw = zeros(size(w));
   
    for j = 1 : numel(wn)
        if (w(j)>=20)
            dw(j)	= pipetting_error(w(j),wn(j));  %the pipetting error is calculated 
        end
    end;
    
    % find and replace zeros in the matrix of volumes
    pz=p;
    f=find(p<1);
    if (f)
          pz(f)=1;
    end
     wz=w;
     f=find(w<1);
    if (f)
         wz(f)=1;
    end
      
   % calculate the new dilution every time, using tmp vars
   for i= 1: length(p(:,1)) %for every dilution step
    gamma_old	= gamma;
    dgamma_old	= dgamma;
    gamma_new	= sum(p(i,:)) / (sum( [ p(i,:) w(i,:)]));
    dgamma_new	= gamma_new * sqrt(	sum( ( [ dp(i,:) dw(i,:) ] ./ [ pz(i,:) wz(i,:) ] ).^2 ) );
    gamma	= gamma_new * gamma_old;
    dgamma	= gamma * sqrt( ( dgamma_new / gamma_new ) ^2 + ( dgamma_old / gamma_old ) ^2 );
   end;

end	% dilution function
  
