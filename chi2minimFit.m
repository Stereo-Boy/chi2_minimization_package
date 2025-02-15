function [bestParams,bestChiSq, yy, SE, stdFits]=chi2minimFit(funVa, pMin, pMax, xVal, yVal, nVal, nbFits,verbose,xx,opts)
% This is a generic function whose goal is to provide a chi square minimization varying the parameters 
% of a function. It loops several times with random starting parameter values to avoid local minima. 
% The best parameters (lowest chi square) are then selected. 
%
% Typically, it is used that way:
% bestParams=chi2minimFit(funVal, pmin, pmax, xVals, yVals, nVals)
%
% Advanced use:
% opts = optimset('display','off');
% [bestParams,bestChiSq,yy,SE,stdFits]=chi2minimFit(funVal, pMin, pMax, xVals, yVals, nVals, nbFits, 'verboseOn', xx,opts)
% ---------------------------------------------------------------------
% Important: this function harnashes lsqnonlin but transforms it in a chi square minimization by calculating a value E 
% in the error function, which gives the chi square value when squared and summed up by lsqnonlin.
% E formula is: E =             ymodel - yVals
%                   ------------------------------------
%                   sqrt(eps+ ymodel.*(1 - ymodel)./nVals)
% with nVals being the number of data to calculate each point in yVals. So nVals is crucial.
% It only works if yVals and ymodel varies between 0 and 1 and it makes more sense when they are probabilities.
% funVal needs to be of the form: prob = function_name(xVals, varying_p) with the different parameter values of the varying ones being in varying_p(1),  varying_p(2)...
% Error is calculated automatically through errorfun so that funVal should NOT calculate error
% Finally, if an xx is provided, the error function funVal should be able to calculate the ymodel (prob) for xx when we provide an yVals input that is 0
% ---------------------------------------------------------------------
% Inputs
% funVal - the function name or function handle of the function to minimze
% pMin - a vector of minimal values for the carying parameters of the function
% pMax - a vector of maximal values
% xVals - a vector of x data to pass to the function
% yVals - a vector of y data to pass to the function
% nVals - a vector of number of trials for each xVals,yVals datapoint
%   nfits - the number of times to loop the minimization to avoid dependence to starting values of the parameters
%   verbose - either verboseON (default) or verboseOFF
%   xx for plot - this is the xx values for which we will calculate the model yy with the best found parameters
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
try
global funVal xVals yVals nVals % this is ugly and not advised to do but lsqnonlin does not allow to pass additionnal fixed parameters that can be then sent to the error function without being varied and we need those
funVal = funVa; xVals = xVal; yVals = yVal; nVals = nVal;
if ~exist('nbFits','var')||isempty(nbFits);nbFits=10; end
if ~exist('verbose','var')||isempty(verbose);verbose='verboseON'; end
if ~exist('xx','var')||isempty(xx);xx=[]; end
if ~exist('opts','var')||isempty(opts);opts=[];
    if strcmp(verbose,'verboseOFF'); opts= optimset('display','off'); end
end
if any(yVals<0 )|| any(yVals>1)
    error('Data are out of 0-1 bounds - manual normalization required.')
end

% add path to functions and define some folder paths, check that they exist
thisPath=fileparts(mfilename('fullpath')); %path to folder where code is located
shared_path = fullfile(thisPath,'shared_functions'); %path to common functions
fun_path = fullfile(thisPath,'fitting_functions'); %path to fitting functions
check_folder(shared_path,1,verbose);  addpath(genpath(shared_path));  % add that path to use these functions
check_folder(fun_path,1,verbose); addpath(genpath(fun_path)); % add that path to use these functions

dispi('Fitting with chi square minimization', verbose)
bestChiSq=inf;
err_fun = @(varying_p)(yVals - feval(funVal,xVals, varying_p))./sqrt(eps+ feval(funVal,xVals, varying_p).*(1 - feval(funVal,xVals, varying_p))./nVals);  %(this formula is from Stanley Klein)
storeParams = nan(nbFits,numel(pMin));
if strcmp(verbose, 'verboseON'); progBar = ProgressBar(nbFits,['Computing ',num2str(nbFits),' fits...']); end
for j=1:nbFits %test different starting parameters and try to get the ones that gives the best chi square
    pStart=rand(1,numel(pMin)).*(pMax-pMin)+pMin;
    [paramsLSQ,chisqLSQ,~,~,~,~,jacobian] = lsqnonlin(err_fun,pStart,pMin,pMax,[],[],[],[],[],opts);
    if chisqLSQ<bestChiSq %select the best chi square parameters here
        bestChiSq=chisqLSQ; bestParams=paramsLSQ; best_j = jacobian;
    end
    storeParams(j,:) = paramsLSQ;
	if strcmp(verbose, 'verboseON'); progBar(j); end
end
stdFits = std(storeParams);

%CALCULATE SE
j=full(best_j);  %something about sparse matrices
covar=inv(j'*j);  %see Numerical Recipes for this formula
SE=sqrt(diag(covar))';%Standard error
                
% calculate the yy model values with best parameters when xx is provided
yy = [];
if ~isempty(xx);    yy = feval(funVal, xx, bestParams); end

dispi('Best parameters : ', bestParams, verbose)
catch err
%% ===== DEBUGING =====%
    disp(err); disp('Error occured: to throw the error, run rethrow(err).')
    disp('Debugging mode')
    keyboard;
end