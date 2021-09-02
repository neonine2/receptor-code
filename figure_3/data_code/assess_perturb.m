function [] = assess_perturb(fname,cellmodel,varargin)
%% testing robustness of optimal strategy to deviation. 
p = inputParser;
addRequired(p,'fname',@isfile);
addRequired(p,'cellmodel',@isnumeric);

% optional arguments
addParameter(p,'saving_data',true,@islogical);
parse(p,fname,cellmodel,varargin{:});
fname = p.Results.fname;
cellmodel = p.Results.cellmodel;
saving_data = p.Results.saving_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% importing optimization result
load(fname,'envmean','receptor_params','conversion_factor',...
                'optr','radlist');
[nloc,m] = size(envmean,[1,2]);
optr = optr(:,:,cellmodel);
envmean = envmean(:,:,cellmodel);
rparams = receptor_params;
rparams.kd = receptor_params.kd*conversion_factor(cellmodel);

% shift and flatten perturbation
shiftparam = linspace(-20,20,21); %should include 0
flattenparam = linspace(1,61,21); %should include 1
perturbparam = combvec(shiftparam,flattenparam)';
nparam = size(perturbparam,1);
perturbMI = zeros(nloc,nparam);
perturbr = zeros(nloc,m,nparam);
parfor ii = 1:nparam
    xtrans = perturbparam(ii,1);
    flatf = perturbparam(ii,2);
    for jj = 1:nloc
        cmean = envmean(jj,:);
        %shift by xtrans and flatten by flatf
        flattened = movmean(optr(jj,:),flatf);
        perturbed = circshift(flattened, xtrans);
        perturbr(jj,:,ii) = perturbed;
        perturbMI(jj,ii) = -totalMI(perturbed,cmean,0,rparams);
    end
    if mod(ii,50) == 0
        disp(ii);
    end
end

cellrad = radlist(cellmodel);
[filepath,name,ext] = fileparts(fname);
if saving_data
    save(strcat(name,"_perturb_",num2str(cellrad),"um"),'perturbMI',...
        'perturbr','perturbparam','shiftparam','flattenparam',...
        'nparam','cellmodel','cellrad')
    disp('Done saving!')
end
end