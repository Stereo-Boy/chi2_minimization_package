function e = errorfun(varying_p,fixed_p)
% varying_p (varying parameters) are parameters that lsqnonlin is trying to fit
% fixed_p is a structure with the following:
%   -  fixed_p.fun - the function that we are trying to estimate
%   -  fixed_p.x fixed_p.y - the x and y to fit
%   -  fixed_p.n - the number of points in each estimate of y
%   -  fixed_p.others - the fixed parameters for the function fun

prob = feval(fixed_p.fun,fixed_p.x, varying_p, fixed_p.others);
e = (fixed_p.y - prob)./sqrt(eps+ prob.*(1 - prob)./fixed_p.n);  %(this line is from Stanley Klein)
