function [ent] = compute_entropy(c)
%Take ligand profile and compute entropy assuming poisson statistics
mfac = 10;
ent = 0;
for ii = 1:length(c)
    mu = c(ii);
    xmax = ceil(mu) * mfac ;
    pdf = poisspdf(0:xmax,mu);
    ent = ent - sum(pdf.*log(pdf),'omitnan');
end

end

