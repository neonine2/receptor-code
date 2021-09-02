function [grad,ex] = approx_grad(cmu,params)
%APPROX_GRAD Summary of this function goes here
%   Detailed explanation goes here
ex = sum(arrayfun(@(x) hillfun(x,params)*poisspdf(x,cmu),...
    0:ceil(cmu)*50));
exlogx = sum(arrayfun(@(x) hillfun(x,params)*log(hillfun(x,params))*poisspdf(x,cmu),...
    0:ceil(cmu)*50));

grad =  exlogx - ex*log(ex);

end

