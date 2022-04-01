clear all;
close all;

cellmodel = 2; % 10um radius

% getting a sequence of environment encountered by a cell being translated
fname = which("soil_var_2_szopt_scheme_h_dynamic_10um.mat");
load(fname,"schemec","scheme_param")
scheme_param.receptornoise = 0.1;
m = size(schemec,2);
epsilon = 0.1/m;
% Optimization parameter for computing optimal placement (in dynamic
% setting)
problem.options = optimoptions('fmincon','Algorithm','interior-point',...
                         'SpecifyObjectiveGradient',true,...
                         'Display','off','CheckGradient',false);
problem.solver = 'fmincon';
problem.Aeq = ones(1,m);
problem.beq = 1;
problem.lb = zeros(1,m);
problem.ub = ones(1,m);
x0 = (-epsilon + (epsilon-(-epsilon)).*rand(1,m) + 1/m)*scheme_param.rtot;
problem.x0 = x0./sum(x0);

ngamma = 11;
gammalist = logspace(-4,-1.5,ngamma);
xoptMAT = zeros(2,m,ngamma); % t
xopt2MAT = zeros(2,m,ngamma); % t = t + delta t
index = [23,53]; % select two representative sample

parfor jj = 1:ngamma
    prob = problem;
    for ii = 1:2
        env = schemec(index(ii),:);
        exlogx = compute_exlogx(env,scheme_param);
        prob.objective = @(x) totalMI(x,env,exlogx,scheme_param);
        [xopt,~] = fmincon(prob);
        xoptMAT(ii,:,jj) = xopt;

        env2 = schemec(index(ii)+1,:);
        exlogx = compute_exlogx(env2,scheme_param);
        prob.objective = @(x) totalMI(x,env2,exlogx,scheme_param,...
            'gamma',gammalist(jj),'lastrvec',xopt);
        [xopt2,~] = fmincon(prob);
        xopt2MAT(ii,:,jj) = xopt2;
    end
end

save("gamma_dynamics.mat")

%% panel S5E
load("gamma_dynamics.mat")
close all
tiledlayout(4,2,'TileSpacing','Compact')
ylimmat = [[0,0.3];[0,0.6]];
for jj = 6:9
    for ii = 1:2
        nexttile
%         title(sprintf('gamma = %0.1e',gammalist(jj)),'fontsize',12)
        yyaxis left
        plot(circshift(schemec(index(ii),:),10),':','color','#D95319','linewidth',1)
        hold on
        plot(circshift(schemec(index(ii)+1,:),10),'-','color','#D95319','linewidth',1)
        ylim(ylimmat(ii,:))
        yyaxis right
        area(circshift(xoptMAT(ii,:,jj),10),'LineStyle',':','linewidth',0.5,...
            'FaceColor',[253 245 157]/255,'Facealpha',0.7)
        area(circshift(xopt2MAT(ii,:,jj),10),'linestyle','-','linewidth',0.5,...
            'FaceColor',[253 245 157]/255,'Facealpha',0.7)
        hold off
        box off
        axis off
        pbaspect([1,1,1])
        ylim([0,0.1])
    end
end
saveas(gca,'panel_S5E.png')

%% 