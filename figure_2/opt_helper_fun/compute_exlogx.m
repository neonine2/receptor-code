function [output] = compute_exlogx(mu,params)

if any(mu < 0)
    error('mean must be non-negative')
end
m = length(mu);
mfac = 10;
xmax = ceil(mu)* mfac;
output = zeros(1,m);

for ii = 1:m
    xm = xmax(ii);
    px = poisspdf(0:xm, mu(ii));
    exlogx = sum(px.*(hillfun(0:xm,params).*log(hillfun(0:xm,params))),'omitnan');
    output(ii) = exlogx;
end

end

