function yprob = logistic_f_chi2(x, varying_p)
% logistic likelihood function used in Serrano-Pedraza et al., 2016
%
% It is a logistic function of probability of correct response as a function of x
% g is guess rate
% lapse is finger error rate 
% delta is the function extent considered for calculation [delta to 1-delta]
% p is the performance level defining threshold (usually 0.75 for a 2AFC)
% x is the x data
% thresholds is the threshold in x units
% pos_slope is the standard deviation in x units

thresholds=varying_p(1);
pos_slope=varying_p(2);
g = 0.5;
delta = 0.01;
p = 0.75;
lapse = 0.05;
b = (2./pos_slope).*log((1 - lapse - g - delta)./delta);
a = (1./b).*log((1-lapse-p)./(p-g));
yprob = g + ((1 - lapse - g)./(1 +exp(-b.*(a + x - thresholds))));

end