clear all;
close all;

cellmodel = 2; % 10um radius

% getting a sequence of environment encountered by a cell being translated
fname = which("soil_var_2_szopt_scheme_h_dynamic_10um.mat");
load(fname,"schemec","scheme_param")
m = size(schemec,2);
epsilon = 0.1/m;
% Optimization parameter for computing optimal placement (in dynamic
% setting)
problem.options = optimoptions('fmincon','Algorithm','sqp',...
                         'SpecifyObjectiveGradient',true,...
                         'Display','off','CheckGradient',false);
problem.solver = 'fmincon';
problem.Aeq = ones(1,m);
problem.beq = 1;
problem.lb = zeros(1,m);
problem.ub = ones(1,m);
x0 = (-epsilon + (epsilon-(-epsilon)).*rand(1,m) + 1/m)*scheme_param.rtot;
problem.x0 = x0./sum(x0);

xoptMAT = zeros(2,m); % t
xopt2MAT = zeros(2,m); % t = t + delta t
index = [23,53]; % select two representative sample
gamma = 0.004;
for ii = 1:length(index)
    env = schemec(index(ii),:);
    exlogx = compute_exlogx(env,scheme_param);
    problem.objective = @(x) totalMI(x,env,exlogx,scheme_param);
    [xopt,~] = fmincon(problem);
    xoptMAT(ii,:) = xopt;

    env2 = schemec(index(ii)+1,:);
    exlogx = compute_exlogx(env2,scheme_param);
    problem.objective = @(x) totalMI(x,env2,exlogx,scheme_param,...
        'gamma',gamma,'lastrvec',xopt);
    [xopt2,~] = fmincon(problem);
    xopt2MAT(ii,:) = xopt2;
end

%% panel 4B
close all
tiledlayout(1,2)
for ii=1:2
    nexttile
    yyaxis left
    plot(circshift(schemec(index(ii),:),10),'k:','linewidth',1)
    hold on
    plot(circshift(schemec(index(ii)+1,:),10),'k-','linewidth',1)
    yyaxis right
    area(circshift(xoptMAT(ii,:),10),'LineStyle',':','linewidth',1,...
        'FaceColor',[253 245 157]/255,'Facealpha',0.7)
    area(circshift(xopt2MAT(ii,:),10),'linestyle','-','linewidth',1,...
        'FaceColor',[253 245 157]/255,'Facealpha',0.7)
    hold off
    box off
    axis off
    pbaspect([1,1,1])
end
saveas(gca,'panel_4B.svg')