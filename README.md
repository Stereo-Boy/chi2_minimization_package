# chi2 minimization package
March 2019 - Adrien Chopin

## Description
This is a generic fitting algorithm whose goal is to provide a chi square minimization technique varying the parameters 
of a given function to fit some data. It loops several times with random starting parameter values to avoid local minima. 
The best parameters (chosen by lowest chi square) are then selected.

## Important things to consider before using
This function harnashes the power of lsqnonlin but transforms it in a chi square minimization by calculating a value E in the error function, which gives the chi square value when squared and summed up by lsqnonlin.
E formula is: E =             ymodel - y
                   ------------------------------------
                   sqrt(eps+ ymodel.*(1 - ymodel)./n)

with n being the number of data to calculate each point in y. So n is crucial. It only works if y and ymodel varies between 0 and 1 and it makes more sense when they are probabilities.

## Usage
### Define your fitting function
You will need to write two functions and store them in the fitting_functions folder:
* the fitting function, used for plotting functionality
* fun, the corresponding chi2 function, whose name ends up with _chi2, used to calculate errors

fun needs to be of the form: prob = your_fitting_function_chi2(x, varying_p, fixed_p) with the different parameter values of the varying ones being in varying_p(1),  varying_p(2)..., and the values of the fixed parameters being in fixed_p(1)...Error is calculated automatically through errorfun so that fun should NOT calculate error. Finally, if an xx is provided to the main call (see below), fun should be able to calculate the ymodel (prob) for xx when we provide an y input.

Only the following functions are currently provided:
* logistic
* probit

### Main call
The main function to call is chi2minimFit.m
Typically, it is used that way:
```
bestParams = chi2minimFit(fun, pmin, pmax, x, y, n)
```
Advanced usage includes:
```
opts= optimset('display','off')
[bestParams,bestChiSq,yy,SE,stdFits]=chi2minimFit(fun, pMin, pMax, x, y, n, other_fixed_p, nbFits, verbose, xx, opts)
```
#### Inputs
* fun - the function name or function handle of the function to minimize (in fitting_functions folder)
* pMin - a vector of minimal values for the each varying parameters of the function
* pMax - a vector of maximal values for the each varying parameters of the function
* x - a vector of x data to pass to the function
* y - a vector of y data to pass to the function
* n - a vector of number of trials for each x,y datapoint
* optional other_fixed_p: fixed parameters values
* optional nfits - the number of times to loop the minimization to avoid dependence to starting values of the parameters
* optional verbose - either verboseON (default) or verboseOFF
* optional xx for plot - this is the xx values for which we will calculate the model yy with the best found parameters
% optional opts - an option structure for lsqnonlin

#### Outputs
* bestParams - the best parameters
* bestChiSq - the value of best chi square (from LSQ 2-norm of squared residuals)
* yy - model values for xx with the best parameters
* SE is the standard error for each parameter using Jacobian estimate
Be careful if you estimate your parameters in a log space: the SE will also be in log space. So to plot it, do: ciplot(10.^(threshold-log10(SE_thr)),10.^(threshold+log10(SE_thr)),xx,'r',0.3);
* stdFits is the std of parameters across repetitions of the different fits - if ~0, you do not need to use a high nbFits
