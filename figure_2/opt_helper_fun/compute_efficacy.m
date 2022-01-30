function [efficacy] = compute_efficacy(optMI,unifMI,varargin)
%COMPUTE_EFFICACY Summary of this function goes here
p = inputParser;
addRequired(p,'unifMI',@isnumeric); % tissue, grad, soil
addRequired(p,'optMI',@isnumeric); 

% optional arguments
addParameter(p,'mode','diff', @isstring);

parse(p,unifMI,optMI,varargin{:});
unifMI = p.Results.unifMI;
optMI = p.Results.optMI;
mode = p.Results.mode;

if isequal(mode,'diff')
    efficacy = (optMI-unifMI)*log2(exp(1));
elseif isequal(mode,'expdiff')
    efficacy = exp(1).^optMI-exp(1).^unifMI;
elseif isequal(mode, 'ratio')
    efficacy = (optMI-unifMI)./unifMI .*100;
elseif isequal(mode, 'unif')
    efficacy = unifMI*log2(exp(1));
end

end

