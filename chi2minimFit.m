function [bestParams,bestChiSq, yy, SE, stdFits]=chi2minimFit(fun, pMin, pMax, x, y, n, other_fixed_p,nbFits,verbose,xx,opts)
% This is a generic function whose goal is to provide a chi square minimization varying the parameters 
% of a function. It loops several times with random starting parameter values to avoid local minima. 
% The best parameters (lowest chi square) are then selected. 
%
% Typically, it is used that way:
% bestParams=chi2minimFit(fun, pmin, pmax, x, y, n)
%
% Advanced use:
% [bestParams,bestChiSq,yy,SE,stdFits]=chi2minimFit(fun, pMin, pMax, x, y, n, other_fixed_p, nbFits, verbose, xx,opts)
% and typical opts is opts= optimset('display','off');
% ---------------------------------------------------------------------
% Important: this function harnashes lsqnonlin but transforms it in a chi square minimization by calculating a value E 
% in the error function, which gives the chi square value when squared and summed up by lsqnonlin.
% E formula is: E =             ymodel - y
%                   ------------------------------------
%                   sqrt(eps+ ymodel.*(1 - ymodel)./n)
% with n being the number of data to calculate each point in y. So n is crucial.
% It only works if y and ymodel varies between 0 and 1 and it makes more sense when they are probabilities.
% Fun needs to be of the form: prob = defineLikelihood_bell_chi2(x, varying_p, fixed_p)
% with the different parameter values of the varying ones being in varying_p(1),  varying_p(2)..., and the values 
% of the fixed parameters being in fixed_p(1)... 
% Error is calculated automatically through errorfun so that fun should NOT calculate error
% Finally, if an xx is provided, the error function Fun should be able
% to calculate the ymodel (prob) for xx when we provide an y input that is 0
% ---------------------------------------------------------------------
% Inputs
% fun - the function name or function handle of the function to minimze
% pMin - a vector of minimal values for the carying parameters of the function
% pMax - a vector of maximal values
% x - a vector of x data to pass to the function
% y - a vector of y data to pass to the function
% n - a vector of number of trials for each x,y datapoint
% other_fixed_p: fixed parameters values
%   nfits - the number of times to loop the minimization to avoid
%       dependence to starting values of the parameters
%   verbose - either verboseON (default) or verboseOFF
%   xx for plot - this is the xx values for which we will calculate the
%       model yy with the best found parameters
%   opts - an option structure for lsqnonlin
% Outputs
% bestParams - the best parameters
% bestChiSq - the value of best chi square (from LSQ 2-norm of squared residuals)
% yy - model values for xx with the best parameters
% SE is the standard error for each parameter using Jacobian estimate
% Be careful if you estimate your parameters in a log space: the SE will also be in log space. 
% So to plot it, do: ciplot(10.^(threshold-log10(SE_thr)),10.^(threshold+log10(SE_thr)),xx,'r',0.3);
% stdFits is the std of parameters across repetitions of the different fits - if ~0, 
% you do not need to use a high nbFits
% ---------------------------------------------------------------------
% Written in March 2019 - Adrien Chopin
% ---------------------------------------------------------------------
% Needs the dispi function
% ---------------------------------------------------------------------

   if ~exist('nbFits','var')||isempty(nbFits);nbFits=10; end
   if ~exist('verbose','var')||isempty(verbose);verbose='verboseON'; end
   if ~exist('xx','var')||isempty(xx);xx=[]; end
   if ~exist('opts','var')||isempty(opts);opts=[]; 
       if strcmp(verbose,'verboseOFF'); opts= optimset('display','off'); end
   end

   if any(y<0 )|| any(y>1)
      error('This function only works for y between 0 and 1 - try normalization to use this function')
   end
   
dispi('Fitting with chi square minimization', verbose)
bestChiSq=inf; fixed_p.fun = fun; fixed_p.x = x; fixed_p.y = y; fixed_p.n = n; fixed_p.others = other_fixed_p;
storeParams = nan(nbFits,numel(pMin));
for j=1:nbFits %test different starting parameters and try to get the ones that gives the best chi square
    pStart=rand(1,numel(pMin)).*(pMax-pMin)+pMin;
    [paramsLSQ,chisqLSQ,~,~,~,~,jacobian] = lsqnonlin('errorfun',pStart,pMin,pMax,opts,fixed_p);
    if chisqLSQ<bestChiSq %select the best chi square parameters here
        bestChiSq=chisqLSQ; bestParams=paramsLSQ; best_j = jacobian;
    end
    storeParams(j,:) = paramsLSQ;
end
stdFits = std(storeParams);

%CALCULATE SE
   j=full(best_j);  %something about sparse matrices
   covar=inv(j'*j);  %see Numerical Recipes for this formula
   SE=sqrt(diag(covar))';%Standard error
                
% calculate the yy model values with best parameters when xx is provided
if ~isempty(xx)
   yy = feval(fun, xx, bestParams, fixed_p.others);
else
   yy=[];
end

dispi('Best parameters : ', bestParams, verbose)