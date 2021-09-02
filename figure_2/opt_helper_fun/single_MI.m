function [MI,efxlog,pxy] = single_MI(rvec,mu,params)

if isnan(rvec)
    rvec = 0;
end

if mu < 0 
    error('ligand mean must be non-negative')
elseif (mu == 0) || (rvec == 0)
    MI = 0;
    efxlog = hillfun(mu,params).*log(hillfun(mu,params));
else
    mfac = 10;
    xmax = ceil(mu) * mfac ;
    if ceil(hillfun(xmax,params)*rvec) > 20
        mfac = 5;
    end
    ymax = ceil(hillfun(xmax,params)*rvec) * mfac;
%     if rvec < 0
%         rvec = 0;
%     end
%     ymax = round(rvec*params.rtot); % binomial channel
    allcomb = combvec(0:xmax,0:ymax)';
    X = allcomb(:,1);
    Y = allcomb(:,2);
    clear allcomb

    fX = hillfun(X,params) * rvec;
    pybarx = poisspdf(Y,fX);
%     prob = hillfun(X,params)/params.rtot; % binomial channel
%     pybarx = binopdf(Y,ymax,prob); % binomial channel
    pxy = pybarx.*poisspdf(X,mu);
    py = sum(reshape(pxy,1+xmax,1+ymax),1)';

    logpybarx = log(pybarx);
    logpybarx(isinf(logpybarx)) = 0;
    logpy = log(py);
    logpy(isinf(logpy)) = 0;

    MI = sum(pxy.*logpybarx) - sum(py.*logpy);

    % for computing gradient
    pxy = reshape(pxy,xmax+1,ymax+1)';
    efxbary = sum(pxy./py .* hillfun(0:xmax,params),2,'omitnan');
    efxlog = sum(py.*(efxbary.*log(efxbary)),'omitnan');
end

end