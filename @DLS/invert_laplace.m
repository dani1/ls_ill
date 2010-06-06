%============================================================================
% FIT CORRELOGRAMMS USING INVERSE LAPLACE TRANSFORM (CONTIN)
%============================================================================
function [ LP_Gamma LP_D LP_Dist ] = invert_laplace ( t, g, q, alpha, cycles )
% The core function is rilt.m by another author. In this function, the data is prepared
% before and after the CONTIN algorithm. Also, CONTIN is performed many times up to a certain
% precision, using low-resolution fits as input for better ones.
%
% NB: the Laplace transform is performed not in terms of decay rates, but of decay times.
 
 [t g]	= filter_correlogramms(t,g);			% eliminate negative correlations

 y	= sqrt(g);					% field correlation function

 smin	= 1/ (1e6 * 100 * q^2);				% lower limit for the s space...
 smax	= 1/ (1e6 * 0.1 * q^2);				% ...and upper limit

 for i = 1 : cycles 					% repeat CONTIN many times with increasing precision
 
  sn	= 10 + 20*(i-1);				% precision for the s space
  s0	= logspace(log10(smin),log10(smax),sn);		% generate initial s space

  if i == 1						% initial guess distribution
   g0 = initial_distr(s0,q);
  else
   g0	= interp1(s1,g,s0,'linear');			% interpolate fitted g to the new space
  end

  ssdthreshold	= max( 1e-8, 1e-5 / 10^i );		% threshold for the fit

  [g,yfit,cfg]	= rilt( t, y, s0, g0, alpha, ssdthreshold );		% perform the fit
  s1	= s0;						% s1 is the old s

 end

 [ LP_Gamma ind ]	= sort(1./s0);			% transform decay times into decay rates...
 LP_D			= 1e-6.*LP_Gamma./q.^2;		% transform decay times into decay rates...
 LP_Dist		= g(ind)';			% ...and reorder the distribution vector

end	% invert_laplace

%============================================================================
% FILTER CORRELOGRAMMS
%============================================================================
function [ t g ] = filter_correlogramms( t, g )

 index1	= g>0;						% retain only positiv gs
 index2	= t< 10;					% maximum time [ms]	

 index	= index1 & index2;				% combine the criteria

 t	= t(index);
 g	= g(index);

end

%============================================================================
% GUESS INITIAL DISTRIBUTION
%============================================================================
function g0 = initial_distr ( s0, q );

     g0 = zeros(length(s0),1);

     D0	= 10;						% A^2/ns

     gamma	= D0 * q^2;				% guessed gamma
     tau	= 1 / gamma;				% guessed tau

     [tmp ind]	= min(s0 - tau);			% find the nearest element of s0 to tau

     varn	= floor(0.1*length(s0));

     for i = 1 : varn
      g0(ind-floor(0.5*varn)+i)	= 1;
     end

     g0	= g0 / sum(g0);					% ...and normalize it

end

%============================================================================
% RILT CONTIN-LIKE ALGORITHM (modified by Fabio Zanini)
%============================================================================
function [g,yfit,cfg] = rilt(t,y,s,g0,alpha,ssdthreshold)
%     rilt   Regularized Inverse Laplace Transform
% 
% [g,yfit,cfg] = rilt(t,y,s,g0,alpha)
% 
% Array g(s) is the Inverse Laplace Transform of the array y(t),
% calculated by a regularized least squares method. This script is
% an emulation of S. Provencher CONTIN program, written in Fortran.
% See http://s-provencher.com/pages/contin.shtml.
% 
% y: data array, a function of the array t (for example: time).
% s: the s space, defined by the user.
% g0: the guess distribution, defined by the user.
% 
% yfit is the Laplace Transform of g and is fitted
% against y at each step of calculation.
% 
% cfg is a structure containing the parameters of the calculation.
% 
% alpha: the regularization parameter (see CONTIN manual,
% by S. Provencher: http://s-provencher.com/pages/contin.shtml,
% Section 3.2.2).
% In this script the second derivative of g is minimized with
% weight alpha^2 ("principle of parsimony"): hor high alpha
% values, the output g distribution will be smooth and
% regular. For alpha==0, the output g distribution will
% be totally 'free'.
% 
% TIP: Start with a g0(s) at low resolution (for example, 10
% points), and obtain a g. Then define an s-space at higher
% resolution and create a corrispondent g0 by interpolating the
% g, for example:
%   g0 = interp1(old_s,g,new_s,'linear');
% Finally repeat the rilt algorithm starting with the new g0.
% 
% REMARKS
% The Inverse Laplace Trasform is a highly ill-posed problem and 
% is therefore intrinsically affected by numerical instability, i.e.
% its solution may not be unique, may not exist or may not depend
% continuously on the data. The g calculated by rilt.m is only
% one of the possible solutions, and this limit is unavoidable.
% Changing the parameters, the solution may change too.
% 
% [...] = rilt(t,y,s,g0,alpha,plot_type,maxsearch,options,shape,constraints,R,w)
% 
% plot_type: defines how to plot the data during the calculation.
% Does not affect the result. Default is 'logarithmic'.
% Available options for plot_type:
%     'logarithmic' - logarithmic x axis
%     'linear' - linear x axis
% 
% maxsearch: Max number of iterations. Default: 1000.
% 
% options: options structure for fminsearch. Default:
%     options = optimset('MaxFunEvals',1e8,'MaxIter',1e8);
%     
% shape: can be 'raise' for raise-up dynamics:
%     y(t) = sum { g(s)*(1-exp(-t/s)) }
% or 'decay' for relaxation dynamics:
%     y(t) = sum { g(s)*exp(-t/s) }
% default is 'decay'.
% 
% constraints: allowed constraints on the g distribution function
% (a cell of multiple constraints can be passed). Available options
% for constraints:
%     'g>0'
%     'zero_at_the_extremes'
% Default is: constraints = {'g>0';'zero_at_the_extremes'};
% 
% R is a matrix determining the form of the Regularizor.
% (see http://s-provencher.com/pages/contin.shtml for details)
% Default settings:
%     R = zeros(length(y)-2,length(y));
% 
% w is the array of fitting weigths for each value of y
% Default settings:
%     w = ones(size(y(:)));
% 
% Leave an empty array to set the default value for the parameters.
% Feel free to contact me and suggest modifications.
% -------------------------------------------------------------
% (c) 2007 Iari-Gabriel Marino, Ph.D.
% University of Parma
% Physics Department
% Viale G.P. Usberti, 7/a
% 43100 Parma - Italy
% e-mail: iarigabriel.marino@fis.unipr.it
% web: http://www.fis.unipr.it/home/marino/
% Tel. +39 (0) 521 906212
% Fax +39 (0) 521 905223
% -------------------------------------------------------------
% 
% Example:
%
% % Target g(s) function:
% g0 = [0 0 10.1625 25.1974 21.8711 1.6377 7.3895 8.736 1.4256 0 0]';
% s0  = logspace(-3,6,length(g0))';
% 
% semilogx(s0,g0); title('Target times distribution - Press any key...'); shg
% pause
% 
% t = linspace(0.01,500,100)';
% [sM,tM] = meshgrid(s0,t);
% A = exp(-tM./sM);
% % Data:
% data = A*g0 + 0.07*rand(size(t));
% 
% alpha = 1e-2;
% 
% % Guess s space and g function
% s = logspace(-3,6,20)';
% g = ones(size(s));
% [g,yfit,cfg] = rilt(t,data,s,g,alpha,'logarithmic',[],[],[],{'g>0'},[],[]);
% 
% subplot(1,2,1)
% plot(t,data,'.',t,yfit,'ro'); title('data and fitting')
% subplot(1,2,2)
% semilogx(s0,g0/max(g0),s,g/max(g),'o-'); title('g-target and g');

%###########################################################

g0 = g0(:);							% must be column
s = s(:);							% must be column
y = y(:);							% must be column
t = t(:);							% must be column

plot_type = @semilogx;						% decide plot type
maxsearch = 100;						% max number of steps
options = optimset('MaxFunEvals',1e8,'MaxIter',1e8);
%constraints = {'g>0','zero_at_the_extremes'};			% constraints	TODO: choose better ones!
constraints = {'g>0'};			% constraints	TODO: choose better ones!
R = zeros(length(g0)-2,length(g0));
w = ones(size(y(:)));
ly = length(y);
oldssd = Inf;

% Prepare the data
[sM,tM] = meshgrid(s,t);
A = exp(-tM./sM);

% Rough normalization of g0, to start with a good guess
g0 = g0*sum(y)/sum(A*g0);

% Plots initialization
fh = gcf; clf(fh);
set(fh,'doublebuffer','on');
s1h = subplot(2,2,1); plot(t,y,'.',t,A*g0); title('Data and fitting curve'); axis tight
s2h = subplot(2,2,2); feval(plot_type,s,g0,'o-'); title('Initial distribution...'); axis tight
s3h = subplot(2,2,3); msdh = plot(0,0); title('Normalized msd');
s4h = subplot(2,2,4); plot(t,abs(y-A*g0)); title('Residuals'); axis tight
msd2plot = 0;
drawnow;

% Main cycle
for k = 1 : maxsearch								% perform the fit at every step, till convergence or limit

 [g,msdfit] = fminsearch(@msd,g0,options,y,A,alpha,R,w,constraints);		% This is the long step: minimize MSD!

 % Re-apply constraints
 if ismember('zero_at_the_extremes',constraints)
     g(1) = 0;
     g(end) = 0;
 end
 if ismember('g>0',constraints)
     g = abs(g);
 end

 g0 = g;									% for the next step
 ssd = sqrt(msdfit/ly);								% Sample Standard Deviation
 ssdStr = num2str(ssd);
 deltassd = oldssd-ssd;								% Difference between "old ssd" and "current ssd"
 oldssd = ssd;

 msd2plot(k) = msdfit/ly;
 plotdata(	s1h,s2h,s3h,s4h,		...
		t,y,A,g,k,			...
		maxsearch,plot_type,s,		...
		msd2plot,msdh,ssdStr,deltassd	);

 if deltassd  < ssdthreshold		break;		end				% Convergence: difference between "old ssd" and "current ssd" < threshold

end

% Saving parameters and results
yfit = A*g; % fitting curve

cfg.t = t;
cfg.y = y;
cfg.yfit = yfit;
cfg.g0 = g0;
cfg.alpha = alpha;
cfg.R = R;
cfg.w = w;
cfg.maxsearch = maxsearch;
cfg.constraints = constraints;
cfg.date = datestr(now,30);

end	% rilt

% ### SUBS #####################################################
function out = msd(g,y,A,alpha,R,w,constraints)
% msd: The mean square deviation; this is the function
% that has to be minimized by fminsearch

% Constraints and any 'a priori' knowledge
if ismember('zero_at_the_extremes',constraints)
    g(1) = 0;
    g(end) = 0;
end
if ismember('g>0',constraints)
    g = abs(g); % must be g(i)>=0 for each i
end
r = diff(diff(g(1:end))); % second derivative of g
yfit = A*g;
% Sum of weighted square residuals
VAR = sum(w.*(y-yfit).^2);
% Regularizor
REG = alpha^2 * sum((r-R*g).^2);
% Output to be minimized
out = VAR+REG;

end	% msd

% ### SUBS #####################################################
function plotdata(s1h,s2h,s3h,s4h,t,y,A,g,k,maxsearch,plot_type,s,msd2plot,msdh,ssdStr,deltassd)
% For the "real-time" plots

axes(s1h);
plot(t,y,'.',t,A*g); title('Data')
xlabel('Time (s)');
axis tight

axes(s2h);
feval(plot_type,s,g,'o-');
title('Relaxation times distribution g(s)');
xlabel('s');
axis tight

axes(s3h);
title(['ssd: ' ssdStr '; \Deltassd: ' num2str(deltassd)])
ylabel('Sample Standard Deviation')
xlabel('Step')
set(msdh,'xdata',1:k,'ydata',msd2plot);
axis tight

axes(s4h);
plot(t,abs(y-A*g)/length(y),'o-');
title('Normalized residuals')
xlabel('Time (s)');
axis tight

set(gcf,'name',['Step:' int2str(k) '/' int2str(maxsearch)])

drawnow

end % plotdata
