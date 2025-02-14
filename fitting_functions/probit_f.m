function y = probit_f(x, m, s, i ,f)
% m - "pivot point" of perf. change, the trial at which you get 50% of the learning
% s - speed of learning - is it a slow learning or a sharp sudden learning at pivot point
% i - initial perf.
% f - final perf.
z = (x-m)./s;
probitf = .5*(erf(z./sqrt(2))+1);%psychometric function going from 0 to 1
y = i + (f-i).*probitf;
