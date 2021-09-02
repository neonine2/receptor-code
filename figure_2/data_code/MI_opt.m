function [rel_eff] = MI_opt(envmodel, fname, radlist, varargin)

p = inputParser;
addRequired(p,'envmodel',@isstring); % tissue, grad, soil
addRequired(p,'fname',@isstring); 
addRequired(p,'radlist',@isnumeric); 

% optional arguments
addParameter(p,'receptor_params',struct('rtot',1000,'kd',40,'receptornoise',0.1), @isstruct); 
addParameter(p,'saving_data',true, @islogical);
addParameter(p,'plot_env',true,@islogical);
addParameter(p,'optgrid',[10,30],@isvector);
addParameter(p,'index','');

parse(p,envmodel,fname,radlist,varargin{:});
envmodel = p.Results.envmodel;
fname = p.Results.fname;
radlist = p.Results.radlist;
receptor_params = p.Results.receptor_params;
saving_data = p.Results.saving_data;
plot_env = p.Results.plot_env;
optgrid = p.Results.optgrid;
index = p.Results.index;

%% Setting parameter values
m = 100; % number of membrane bins

%% setting up environment
disp(strcat('Setting up environment...',fname))
if isequal(envmodel,'tissue') || isequal(envmodel,'grad')
    load(fname,'cbound','csol','xmin','xmax','ymin','ymax');
    ctot = cbound(5:end,:) + csol(5:end,:); % move away from source (left) boundary
    posmat = combvec(linspace(1,xmax-xmin,size(ctot,1)),...
                        linspace(1,ymax-ymin,size(ctot,2)))';
    fenv = scatteredInterpolant(posmat,ctot(:),'natural','linear');
    fconc = @(x,y) fenv([x,y]); % concentration function
    %fitting gradient
    [xmin,xmax,ymin,ymax] = deal(20,250,20,800); %um
    posmat = combvec(1:xmax,1:ymax)';
    env = reshape(fconc(posmat(:,1),posmat(:,2)),xmax,ymax)';
    f = fit((1:xmax)',mean(env)','exp1'); % concentration function
    %tissue environment
    if plot_env && isequal(envmodel,'tissue')   
        imagesc(env) % shows a piece of the entire environment
        title(['mean conc = ',num2str(mean(env,'all'))])
        colorbar()
        figure(2)
        plot(1:xmax,mean(env))
    elseif plot_env && isequal(envmodel,'grad')
        plot(1:xmax,f(1:xmax),1:xmax,mean(env));
        title(['mean conc = ',num2str(mean(env,'all'))])
    end
    pbaspect([1,3,1])
    pause(0.01)
elseif isequal(envmodel,'soil')
    load(fname);
%     clear fconc
    [xmin,xmax,ymin,ymax] = deal(xwin(1)+50,xwin(2)-50,...
                                 ywin(1)+50,ywin(2)-50); %um
    cvar = 9^2; % variance of profile from a single bacterium
    avgconc = 0.6; %nM
    % assume soil density of 1.3 g cm−3 (Raynaud 2014); AHL concentration
    % between 1ug per kg soil, this gives ~1nM, we assume 0.6nM
    if ~exist('fenv') % slow to run
        disp('generating fconc')
        avgapprox = size(soil_LGCP,1)/((xmax-xmin)*(xmax-xmin));
        fconc = @(x,y) sum(mvnpdf([x,y]-soil_LGCP,[0,0],cvar*eye(2)))/avgapprox*avgconc;
        posmat = combvec(xmin:xmax,ymin:ymax)';
        env = arrayfun(fconc,posmat(:,1),posmat(:,2));
        fenv = scatteredInterpolant(posmat,env,'natural','linear');
        fconc = @(x,y) fenv([x,y]); % concentration function
        save(fname,'fconc','fenv','cvar','avgconc','-append')
    end
    if plot_env
        posmat = combvec(xmin+101:xmin+200,ymin+101:ymin+200)';
        env = reshape(fconc(posmat(:,1),posmat(:,2)),100,100)';
        imagesc(env) % shows a piece of the entire environment
        title(['mean conc = ',num2str(mean(env,'all'))])
        colorbar()
        pbaspect([1,1,1])
        pause(0.01)
    end
end

%% setting up cell and cell locations
% setting number of sampled cell locations
lx = optgrid(1); %number of locations in x direction
ly = optgrid(2); %number of locations in y direction
nloc = lx*ly; %total number of sample locations

% setting cell sizes
nrad = length(radlist); % number of cell radius
angle = linspace(0,2*pi*(1-1/m),m)';
cellsurf = reshape(radlist,[1,1,nrad]).*[cos(angle),sin(angle)];

% determining lattice of cell location and conversion factor
centerlist = zeros(nloc,2,nrad);
conversion_factor = zeros(nrad,1);
for jj = 1:nrad
    cellcoord = cellsurf(:,:,jj);
    r = ceil(max(abs(cellcoord),[],'all'));
    x = linspace(xmin+r+1,xmax-r-1,lx); % ensure cell boundary within domain
    y = linspace(ymin+r+1,ymax-r-1,ly);
    centerlist(:,:,jj) = combvec(x,y)';
    conversion_factor(jj) = conc2count(radlist(jj),m);
end

%% Optimization procedure
%Optimization parameters
options = optimoptions('fmincon','Algorithm','sqp',...
                       'SpecifyObjectiveGradient',true,...
                       'Display','off',...
                       'CheckGradients',false,...
                       'OptimalityTolerance',1e-10,...
                       'StepTolerance',1e-12); 
                   %for binomial channel, set larger step size
                   %set FiniteDifferenceStepSize to 1e-3
problem.options = options;
problem.solver = 'fmincon';
problem.Aeq = ones(1,m);
problem.beq = 1;
problem.lb = zeros(1,m);
problem.ub = ones(1,m);
epsilon = 0.02./m;

%pre-allocate space
unifMI = zeros(nloc,nrad);
optMI = zeros(nloc,nrad);
optr = zeros(nloc,m,nrad);
envmean = zeros(nloc,m,nrad);

disp('Running optimization...')
for jj = 1:nrad
    rparams = receptor_params;
    cfactor = conversion_factor(jj);
    rparams.kd = receptor_params.kd * cfactor;
    
    % looping over all cell locations
    parfor ii = 1:nloc
        problem_opt = problem;
        center = centerlist(ii,:,jj);
        poscoord = center + cellsurf(:,:,jj);
        %mean ligand counts (cmean)
        if isequal(envmodel,'grad')
            xy = combvec(1:100,...
                            linspace(center(:,2)-10,center(:,2)+10,10))';
            env = reshape(fconc(xy(:,1),xy(:,2)),100,10)';
            f = fit((1:100)',mean(env)','exp1');
            cmean = f(poscoord(:,1))'*cfactor;
        else
            cmean = arrayfun(fconc,poscoord(:,1),poscoord(:,2))*cfactor; 
        end
        if min(cmean)<0
            warning(['location ii = ',num2str(ii),' ',num2str(jj),' has negative mean'])
            continue
        end
        envmean(ii,:,jj) = cmean;
        exlogx = compute_exlogx(cmean,rparams); %for analytical gradient
        problem_opt.objective = @(x) totalMI(x,cmean,exlogx,rparams);
        % random initialization
        startpt = -epsilon + (epsilon-(-epsilon)).*rand(1,m) + 1/m;
        problem_opt.x0 = startpt./sum(startpt);
        % obtaining optimal solution
        [x,~] = fmincon(problem_opt);
        x(x<0) = 0;
        x = x/sum(x);
        
        % storing solution
        unifMI(ii,jj) = -totalMI(ones(1,m)./m,cmean,exlogx,rparams);
        optMI(ii,jj) = -totalMI(x,cmean,exlogx,rparams);
        optr(ii,:,jj) = x; %exact optimal
        if mod(ii,100)==0
            disp(ii);
        end
    end
    % display result
    disp((sum(optMI) - sum(unifMI))./sum(unifMI) .* 100);
end 
disp('Finished')
rel_eff = (sum(optMI) - sum(unifMI))./sum(unifMI) .* 100;
%% saving data
if saving_data
    if nrad > 1
        sz = 'sz';
    else
        sz = '';
    end
    if isequal(envmodel,'grad')
        filename = strcat(fname,'_grad_',sz,'opt',num2str(index));
    else
        filename = strcat(fname,'_',sz,'opt',num2str(index));
    end
    save(filename,'envmean','optMI','unifMI','optr','problem',...
        'conversion_factor',...
            'receptor_params','cellsurf','centerlist',...
            'radlist','fconc','fenv','rel_eff','-v7.3');
    disp('Done saving!')
end

end