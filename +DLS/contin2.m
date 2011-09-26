%============================================================================
% FIT CORRELOGRAMMS USING INVERSE LAPLACE TRANSFORM (CONTIN)
%============================================================================
function [ s g A ] = contin2 ( t, y, dy, smin, smax, m, alpha, cycles )
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

  plotflag = 'no';
  [g,yfit,cfg]	= rilt( t, y, dy, s0, g0, alpha, ssdthreshold, plotflag );		% perform the fit using the Matlab function

  s1	= s0;						% s1 is the old s

 end							% the variables at this point are already s and g

 s = s0;
 A = cfg.A;

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

     [tmp ind]	= min(s0 - tau);			% find the nearest element of s0 to tau
     g0(ind)	= 1;

     go(floor(length(g0)/2)) = 1;
     g0	= g0 / sum(g0);					% ...and normalize it

end

%============================================================================
% RILT CONTIN-LIKE ALGORITHM (modified by Fabio Zanini)
%============================================================================
function [g,yfit,cfg] = rilt(t,y,dy,s,g0,alpha,ssdthreshold, plotflag)
%     rilt   Regularized Inverse Laplace Transform
% 
% [g,yfit,cfg] = rilt(t,y,s,g0,alpha, ssdthreshold)
%
% The fit algorithm is constrained by
% 1. g >= 0
% 2. g(1) = 0 and g(end) = 0

 if nargin < 8
  plotflag = 'no'
 end
 
 g0 = g0(:);							% must be column
 s = s(:);							% must be column
 t = t(:);							% must be column
 y = y(:);							% must be column
 dy = dy(:);							% must be column
 w = 1 ./ dy.^2;
 
 plot_type = @semilogx;						% decide plot type
 maxsearch = 100;						% max number of steps
 R = zeros(length(g0)-2,length(g0));
 ly = length(y);
 oldssd = Inf;
 
 % Prepare the data
 [sM,tM] = meshgrid(s,t);
 A = exp(-tM./sM);
 
 % Rough normalization of g0, to start with a good guess
 g0 = g0*sum(y)/sum(A*g0);
 
 % Plots initialization
 if strcmp(plotflag,'yes')
  fh = gcf; clf(fh);
  set(fh,'doublebuffer','on');
  s1h = subplot(2,2,1); semilogx(t,y,'.',t,A*g0); title('Data and fitting curve'); axis tight;
  s2h = subplot(2,2,2); feval(plot_type,s,g0,'o-'); title('Initial distribution...'); axis tight
  s3h = subplot(2,2,3); msdh = plot(0,0); title('Normalized msd');
  s4h = subplot(2,2,4); plot(t,abs(y-A*g0)); title('Residuals'); axis tight
  msd2plot = 0;
  drawnow;
 end
 
 % Main cycle
 for k = 1 : maxsearch								% perform the fit at every step, till convergence or limit
 
 % options = optimset('MaxFunEvals',1e8,'MaxIter',1e8,'Display','off');
 % [g,msdfit] = fminsearch(@msd,g0,options);	
 
  MA = (-1)* eye(length(g0));
  LB = zeros(1,length(g0));	
  UB = Inf(1,length(g0));
  extr0 = zeros(length(g0)); extr0(1) = 1; extr0(end)=1;
  options = optimset('Algorithm','active-set','Display','off');
  [ g, msdfit ] = fmincon(@msd,g0,MA,LB,extr0,LB,[],[],[],options);
 
  g0 = g;									% for the next step
  ssd = sqrt(msdfit/ly);								% Sample Standard Deviation
  ssdStr = num2str(ssd);
  deltassd = oldssd-ssd;								% Difference between "old ssd" and "current ssd"
  oldssd = ssd;
 
  if strcmp(plotflag,'yes')
   msd2plot(k) = msdfit/ly;
   plotdata(	s1h,s2h,s3h,s4h,		...
  		t,y,A,g,k,			...
  		maxsearch,plot_type,s,		...
  		msd2plot,msdh,ssdStr,deltassd	);
  end
 
  % Convergence: difference between "old ssd" and "current ssd" < threshold
  if deltassd  < ssdthreshold
   break;
  end
 
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
cfg.A = A;
cfg.maxsearch = maxsearch;
cfg.date = datestr(now,30);

 % ### NESTED FUNCTION TO BE MINIMIZED #############################
 function out = msd(g)
 % msd: The mean square deviation; this is the function
 % that has to be minimized by fminsearch
  
 % Regularizator				% The REG is the second derivative??
 %r = diff(diff(g(1:end)));			% second derivative of g
 %REG = alpha^2 * sum((r-R*g).^2);
 r = diff(g(1:end));				% first derivative of g
 REG = alpha^2 * sum(r.^2);
 %REG = alpha^2 * sum( g.^2 );			% squared value of g
 
 % Sum of weighted square residuals
 yfit = A*g;
 VAR = sum(w.*(y-yfit).^2);

 out = VAR+REG;					% Output to be minimized
 
 end	% msd

end	% rilt

% ### SUBS #####################################################
function plotdata(s1h,s2h,s3h,s4h,t,y,A,g,k,maxsearch,plot_type,s,msd2plot,msdh,ssdStr,deltassd)
% For the "real-time" plots

axes(s1h);
semilogx(t,y,'.',t,A*g); title('Data')
xlabel('Time (ms)');
axis tight

axes(s2h);
feval(plot_type,s,g,'o-');
title('Relaxation times distribution g(ms)');
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
