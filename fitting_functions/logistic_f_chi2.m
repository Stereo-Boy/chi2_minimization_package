function prob = logistic_f_chi2(x, varying_p, fixed_p)
% likelihood function used in Serrano-Pedraza et al., 2016
%
% It is a logistic function of probability of correct response as a
% function of log-disparity with a linear return to chance after two
% threshold values
% g is guess rate
% lapse is finger error rate 
% pos_slope is the standard deviation
% delta is the function extent considered for calculation [delta to 1-delta]
% p is the performance level defining threshold (usually 0.75 for a 2AFC)
% disparities is the disparity in log arcsec
% thresholds is the threshold in log arcsec

thresholds=varying_p(1);
pos_slope=varying_p(2);

g = fixed_p(1);
delta = fixed_p(2);
p = fixed_p(3);
lapse = fixed_p(4);

prob = logistic_f(g, lapse, pos_slope, delta, p, x, thresholds);


end