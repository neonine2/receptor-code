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
% saveas(gca,'panel_4B.svg')

%% new panel 4B figure
close all
index = [23,53,24,54];
opt = [xoptMAT;xopt2MAT];
tiledlayout(2,2,'TileSpacing','Compact');
for ii=1:4
    nexttile
    yyaxis left
    plot(circshift(schemec(index(ii),:),10),'linewidth',1)
    yyaxis right
    area(circshift(opt(ii,:),10),'linewidth',0.8,...
        'FaceColor',[253 245 157]/255,'Facealpha',0.7)
    hold off
    box off
    axis off
    pbaspect([2,1,1])
end

%% NEW PANEL 4B code
% clear all;
% close all;
% 
% cellmodel = 2; % 10um radius
% load('scheme_parameter','scheme_param');
% index = [23,53]; % select two representative sample
% % getting a sequence of environment encountered by a cell being translated
% feedback_scheme_efficiency(which("soil_var_2_szopt.mat"),scheme_param,...
%                       "h","dynamic",cellmodel,"selectind",[index,index+1]);
% load("soil_var_2_szopt_scheme_h_dynamic_10um_selectind.mat",...
%                     "schemec","schemecellp","scheme_param");
% 
% scheme_param.receptornoise = 0.1;
% m = size(schemec,2);
% epsilon = 0.1/m;
% % Optimization parameter for computing optimal placement (in dynamic
% % setting)
% problem.options = optimoptions('fmincon','Algorithm','interior-point',...
%                          'SpecifyObjectiveGradient',true,...
%                          'Display','off','CheckGradient',false);
% problem.solver = 'fmincon';
% problem.Aeq = ones(1,m);
% problem.beq = 1;
% problem.lb = zeros(1,m);
% problem.ub = ones(1,m);
% x0 = (-epsilon + (epsilon-(-epsilon)).*rand(1,m) + 1/m)*scheme_param.rtot;
% problem.x0 = x0./sum(x0);
% 
% xoptMAT = zeros(2,m); % t
% xopt2MAT = zeros(2,m); % t = t + delta t
% 
% gamma = 0.004;
% fcount = scheme_param.fcount;
% cellboundary = scheme_param.cellboundary;
% 
% tiledlayout(5,2)
% for ii = 1:10
%     nexttile
%     coord = cellcoord+[0,10]*ii;
%     plot(fcount(coord(:,1),coord(:,2)))
%     coord = cellcoord+[0,10]*(ii+1);
%     hold on
%     plot(fcount(coord(:,1),coord(:,2)))
%     hold off
% %     pause(1)
% end
% 
% for ii = 1:length(index)
%     cellcoord = schemecellp(ii,:) + cellboundary;
%     env = fcount(schemecellp(ii,1),schemecellp(ii,2));
%     exlogx = compute_exlogx(env,scheme_param);
%     problem.objective = @(x) totalMI(x,env,exlogx,scheme_param);
%     [xopt,~] = fmincon(problem);
%     xoptMAT(ii,:) = xopt;
% 
%     env2 = schemec(ii+2,:);
%     exlogx = compute_exlogx(env2,scheme_param);
%     problem.objective = @(x) totalMI(x,env2,exlogx,scheme_param,...
%         'gamma',gamma,'lastrvec',xopt);
%     [xopt2,~] = fmincon(problem);
%     xopt2MAT(ii,:) = xopt2;
% end
