function [output,grad,summarystat] = totalMI(rvec,xmean,exlogx,params,varargin)

p = inputParser;
isvecorempty = @(x) isvector(x) || isempty(x);
addRequired(p,'rvec',@isvector);
addRequired(p,'xmean',@isvector);
addRequired(p,'exlogx',isvecorempty);
addRequired(p,'params',@isstruct);

% optional arguments
addParameter(p,'lastrvec',1,@isvector);
addParameter(p,'gamma',0,@isnumeric);
addParameter(p,'poiscor',0,@isnumeric);
addParameter(p,'gausscor',0,@isnumeric);

parse(p,rvec,xmean,exlogx,params,varargin{:});
rvec = p.Results.rvec;
xmean = p.Results.xmean;
exlogx = p.Results.exlogx;
params = p.Results.params;
lastrvec = p.Results.lastrvec;
gamma = p.Results.gamma;
poiscor = p.Results.poiscor;
gausscor = p.Results.gausscor;

%compute total mutual information
m = length(xmean);
MI = 0;
exloglist = zeros(1,m);

if poiscor ~= 0
    [MI,exloglist,summarystat] = pairwise_MI(rvec,xmean,poiscor,gausscor,params);
else
    for ii = 1:m
        [MIscalar, efxlogbary] = single_MI(rvec(ii),xmean(ii),params);
        MI = MI + MIscalar;
        exloglist(ii) = efxlogbary;
    end
end  

if gamma > 0 % including wasserstein distance
    [cost, costgrad] = W1dist(rvec,lastrvec);
    MI = MI - gamma*cost;
end
output = -MI;

% gradient only works for independent env
if nargout > 1
    grad = exlogx - exloglist;
    grad(xmean==0) = 0;
    if gamma > 0
        grad = grad(:) - gamma*costgrad(:); % including wasserstein cost
    end
    grad = -grad;
end

if nargout > 2 && poiscor == 0 % gradient required
    warning('no summary statistics available for independent environment')
end

end

