%============================================================================
% FIT CORRELOGRAMMS USING INVERSE LAPLACE TRANSFORM (CONTIN)
%============================================================================
function [ s g ] = contin2 ( t, y, dy, smin, smax, m, alpha, cycles )
% The core function is rilt.m by another author. In this function, the data is prepared
% before and after the CONTIN algorithm. Also, CONTIN is performed many times up to a certain
% precision, using low-resolution fits as input for better ones.
%
% NB: the Laplace transform is performed not in terms of decay rates, but of decay times

 for i = 1 : cycles 					% repeat CONTIN many times with increasing precision
 
  sn	= m + 20*(i-1);					% precision for the s space
  s0	= logspace(log10(smin),log10(smax),sn);		% generate initial s space

  if i == 1						% initial guess distribution
   g0 = initial_distr(s0);
  else
   g0	= interp1(s1,g,s0,'linear');			% interpolate fitted g to the new space
  end

  ssdthreshold	= 1e-13;		% threshold for the fit

  [g,yfit,cfg]	= rilt( t, y, dy, s0, g0, alpha, ssdthreshold );		% perform the fit using the Matlab function

  s1	= s0;						% s1 is the old s

 end							% the variables at this point are already s and g

 s = s0;

end	% invert_laplace

%============================================================================
% GUESS INITIAL DISTRIBUTION
%============================================================================
function g0 = initial_distr ( s0 );

     g0([1:length(s0)]) = 0.01;

     D0	= 1;						% A^2/ns
     q0	= 1e-3;						% A^{-1}

     gamma	= 1e6 * D0 * q0^2;			% guessed gamma
     tau	= 1 / gamma;				% guessed tau

     [tmp ind]	= min(abs(s0 - tau));			% find the nearest element of s0 to tau
     g0(ind)	= 1;

     g0(floor(length(g0)/2)) = 1;
     g0	= g0 / sum(g0);					% ...and normalize it

end

%============================================================================
% RILT CONTIN-LIKE ALGORITHM (modified by Fabio Zanini)
%============================================================================
function [g,yfit,cfg] = rilt(t,y,dy,s,g0,alpha,ssdthreshold)
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
%###########################################################
g0 = g0(:);							% must be column
s = s(:);							% must be column
t = t(:);							% must be column
y = y(:);							% must be column
dy = dy(:);							% must be column
w = 1 ./ dy.^2;

plot_type = @semilogx;						% decide plot type
maxsearch = 100;						% max number of steps
options = optimset('MaxFunEvals',1e8,'MaxIter',1e8,'Display','off');
R = zeros(length(g0)-2,length(g0));
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

 [g,msdfit] = fminsearch(@msd,g0,options);					% This is the long step: minimize MSD!

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
cfg.date = datestr(now,30);

 % ### NESTED FUNCTION TO BE MINIMIZED #############################
 function out = msd(g)
 % msd: The mean square deviation; this is the function
 % that has to be minimized by fminsearch
 
 g(1) = 0;	g(end) = 0;			% G is zero at the extremes
 g = abs(g);					% G is positive
 
 % Regularizator				% The REG is the second derivative??
 r = diff(diff(g(1:end)));			% second derivative of g
 REG = alpha^2 * sum((r-R*g).^2);
 %r = diff(g(1:end));				% first derivative of g
 %REG = alpha^2 * sum(r.^2);
 %REG = alpha^2 * g * g';			% squared value of g
 
 % Sum of weighted square residuals
 yfit = A*g;
 VAR = sum(w.*(y-yfit).^2);

 out = VAR+REG;					% Output to be minimized
 
 end	% msd

 waitforbuttonpress;
 close;

end	% rilt

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
