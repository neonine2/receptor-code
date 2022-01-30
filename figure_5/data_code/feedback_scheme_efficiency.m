function [relative_eff] = feedback_scheme_efficiency(fname,scheme_param,...
                                    testing_param,mode,cellmodel,varargin)
p = inputParser;
addRequired(p,'fname',@isfile); % environment file
addRequired(p,'scheme_param',@isstruct); % environment file
addRequired(p,'testing_param',@isstring); % h or koff or both (dual)
addRequired(p,'mode',@isstring); % static or dynamic
addRequired(p,'cellmodel',@isnumeric); % sets the cell size

% optional arguments
addParameter(p,'saving_data',true,@islogical);
addParameter(p,'selectind',-1,@isvector);

parse(p,fname,scheme_param,testing_param,mode,cellmodel,varargin{:});
fname = p.Results.fname;
scheme_param = p.Results.scheme_param;
testing_param = p.Results.testing_param;
mode = p.Results.mode;
cellmodel = p.Results.cellmodel;
saving_data = p.Results.saving_data;
selectind = p.Results.selectind;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decoder_method = "yaxis1d"; % if static, this option is not needed

% load env
load(fname,'centerlist','envmean','unifMI','optMI','fconc','optr',...
    'conversion_factor','receptor_params','radlist');
envmean = envmean(:,:,cellmodel);
centerlist = centerlist(:,:,cellmodel);
unifMI = unifMI(:,cellmodel);
optMI = optMI(:,cellmodel);
optr = optr(:,:,cellmodel);

% set scheme parameter
scheme_param.mean_cell_radius = radlist(cellmodel);
fcount = @(x,y) fconc(x,y)*conversion_factor(cellmodel);
scheme_param.fcount = fcount;
nminute = 10; % 10 minutes
scheme_param.T = 60*nminute; 
scheme_param.kd = receptor_params.kd*conversion_factor(cellmodel);
scheme_param.rtot = receptor_params.rtot;
scheme_param.receptornoise = receptor_params.receptornoise;
scheme_param.h = 0.004;
scheme_param.koff = 0.18;
time_size = floor(scheme_param.T/scheme_param.dt);

% establish cell shape
m = scheme_param.N;
phi = linspace(pi,-pi+(2*pi)/m,m)';
phi = circshift(flipud(phi),m/2+1);
phi(1) = 0; % correct for minor numerical inaccuracy
s = scheme_param.mean_cell_radius;
scheme_param.cellboundary = s*[cos(phi),sin(phi)];
scheme_param.dx = ellipsePerimeter(s,s)/m;

% parameters to search through
nloc = size(centerlist,1);
if isequal(mode,'static')
    posind = [1:10:nloc,(1:10:nloc)+1];
    if selectind > 0
        posind = posind(selectind);
    end
elseif isequal(mode,'dynamic')
    % make sure cell is far enough from boundary
    posind = 1:10:nloc;
    if nloc > 1000
        posind = 1:40:nloc;
    end
    posind = [posind(4:end),posind(4:end)+1];
    if selectind > 0
        posind = posind(selectind);
    end
end
nparam = 30;
if isequal(testing_param,"h")
    paramlist = logspace(-3,-2,nparam);
elseif isequal(testing_param,"koff")
    paramlist = logspace(-2,-0.4,nparam);
elseif isequal(testing_param,"dual")
    hlist = logspace(-3,-1,16);
    kofflist = logspace(-2,0,16);
    paramlist = combvec(hlist,kofflist)';
    nparam = size(paramlist,1);
end
npos = length(posind);
schemer = zeros(npos,scheme_param.N,nparam);
schemeMI = zeros(npos,nparam);
schemec = zeros(npos,scheme_param.N);
schemecellp = zeros(npos,2);
recstat = zeros(npos,4,nparam);

%% running feedback scheme for different parameter and computing MI 

for ii = 1:nparam
    % set rate parameters
    if isequal(testing_param,"h")
        scheme_param.h = paramlist(ii);
    elseif isequal(testing_param,"koff")
        scheme_param.koff = paramlist(ii);
    elseif isequal(testing_param,"dual")
        scheme_param.h = paramlist(ii,1);
        scheme_param.koff = paramlist(ii,2);
    end
    
    parfor jj = 1:npos
        sparam = scheme_param;
        % set parameter
        sparam.cellp = centerlist(posind(jj),:);
        if isequal(mode,'dynamic') % move back cell
            sparam.cellp = sparam.cellp-[0,sparam.stepsz*2*nminute-1];
        end
        results = cell_move(sparam,"mode",mode,...
                 "decoder_method",decoder_method,...
                 "receptor","feedback");
        % using final receptor placement
        rvec = results.f(end,:)/sum(results.f(end,:)); 
        cmean = results.env(end,:); 
        schemec(jj,:) = cmean;
        schemecellp(jj,:) = results.cellp(end,:);
        schemer(jj,:,ii) = rvec;
        recstat(jj,:,ii) = results.stat;
        % computing MI of scheme-induced receptor placement
        exlogx = compute_exlogx(cmean,sparam);
        schemeMI(jj,ii) = -totalMI(rvec,cmean,exlogx,sparam);
        
%         % observe receptor dynamic
%         for kk = 1:time_size
%             yyaxis left
%             plot(results.f(kk,:)/sum(results.f(kk,:)));
%             yyaxis right
%             plot(results.env(kk,:))
%             pause(0.000001)
%         end
%         % checking ligand profiles match up
%         yyaxis right
%         plot(envmean(posind(jj),:),'Marker','*');
%         hold on;
%         plot(cmean,'k','linewidth',2);hold off;
%         yyaxis left;
%         plot(optr(posind(jj),:),'Marker','*')
%         hold on;
%         plot(rvec,'r','linewidth',2);hold off;
%         pause(1)
    end
    if mod(ii,5) == 0
        disp(ii)
    end
end
disp(squeeze(mean(recstat)));
meankfb = squeeze(mean(recstat(:,1,:)));

% compute scheme-induced distribution efficiency
unifMI = unifMI(posind);
optMI = optMI(posind);
scheme_eff = (sum(schemeMI)-sum(unifMI))./sum(unifMI) * 100;
opt_eff = (sum(optMI)-sum(unifMI))./sum(unifMI) * 100;
relative_eff = scheme_eff/opt_eff;
disp([meankfb,relative_eff']);

cellrad = radlist(cellmodel);
if saving_data
    [~,name,~] = fileparts(fname);
    schemename = strcat("_scheme_",testing_param,"_",mode,"_",num2str(cellrad),"um");
    if selectind > 0
        schemename = strcat(schemename,"_selectind");
    end
    save(strcat(name,schemename),'scheme_param','posind','schemer',...
        'recstat','cellmodel','schemeMI','cellrad',...
        'scheme_eff','opt_eff','relative_eff','schemec','paramlist',...
        'relative_eff','schemecellp');
end

end