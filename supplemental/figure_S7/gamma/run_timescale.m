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
                         'Display','off','CheckGradient',false,...
                         'OptimalityTolerance',1e-12,...
                         'StepTolerance',1e-20);
problem.solver = 'fmincon';
problem.Aeq = ones(1,m);
problem.beq = 1;
problem.lb = zeros(1,m);
problem.ub = ones(1,m);
x0 = (-epsilon + (epsilon-(-epsilon)).*rand(1,m) + 1/m)*scheme_param.rtot;
problem.x0 = x0./sum(x0);

stepsz = 2;
tvec = stepsz:stepsz:30;
Ttot = length(tvec)+1;
ngamma = 11;
gammalist = logspace(-4,-1.5,ngamma);
optMAT = zeros(Ttot,m,ngamma); % t
envMAT = zeros(Ttot,m,ngamma);

fname = which("tissue_300by900_grad_szopt.mat");
load(fname,"envmean")
env = circshift((envmean(1,:,2)'-0.2)*2.5,ceil(Ttot/2));
exlogx = compute_exlogx(env,scheme_param);
problem.objective = @(x) totalMI(x,env,exlogx,scheme_param);
[opt0,~] = fmincon(problem);

parfor jj = 1:ngamma
    xopt = opt0;
    prob = problem;
    gamma = gammalist(jj);
    for ii = 2:Ttot
        cvec = circshift(env,-tvec(ii-1));
        exlogx = compute_exlogx(cvec,scheme_param);
        prob.objective = @(x) totalMI(x,cvec,exlogx,scheme_param,...
                                     'gamma',gamma,'lastrvec',xopt);
        [xopt,~] = fmincon(prob);
        optMAT(ii,:,jj) = xopt;
        envMAT(ii,:,jj) = cvec;
    end
end

for jj = 1:ngamma
    optMAT(1,:,jj) = opt0;
    envMAT(1,:,jj) = env;
end

% compute speed for each gamma
phi = linspace(pi,-pi+(2*pi)/m,m)';
phi = circshift(flipud(phi),m/2+1);
phi(1) = 0; % correct for minor numerical inaccuracy
cellspeed = 2;
%shifted ligand peak by 1 bin per timestep
celldistmoved = stepsz*(Ttot-1) * (2*pi*10/m);
timemoved =  celldistmoved / cellspeed; 
speedMAT = zeros(ngamma,1);
for jj = 1:ngamma
%     [~,I] = max(optMAT([1,Ttot],:,jj),[],2);
    rvec = optMAT([1,Ttot],:,jj);
    z1 = sum(cos(phi).* rvec');
    z2 = sum(sin(phi).* rvec');
    I = atan2(z2,z1);
    angle = min(2*pi-abs(I(2)-I(1)), abs(I(2)-I(1)));
    speedMAT(jj) = (angle/(2*pi) * 2*pi*10)/timemoved;
end
save("gamma_timescale.mat")

%plotting
load("gamma_timescale.mat")
tiledlayout(3,2)
nexttile
plot(env,'--','color','#D95319','linewidth',1)
hold on
plot(envMAT(end,:,1),'-','color','#D95319','linewidth',1)
hold off
title('Ligand profile','FontWeight','normal')
ylabel('ligand level')
xlabel('cell surface')
legend('1 min', '10 min','Location','southeast')
pbaspect([3.5,1,1])
set(gca,'linewidth',1,'fontsize',12)
tile = [3,5];
gammaind = [4,9];
for ii = 1:2
    nexttile(tile(ii))
%     plot(optMAT(1,:,gammaind(ii))'*100,'k--','linewidth',1)
%     plot(optMAT(Ttot,:,gammaind(ii))'*100,'k-','linewidth',1)
    area(optMAT(1,:,gammaind(ii))'*100,'LineStyle',':','linewidth',1,...
            'FaceColor',[253 245 157]/255,'Facealpha',0.7)
    hold on
    area(optMAT(Ttot,:,gammaind(ii))'*100,'linestyle','-','linewidth',1,...
            'FaceColor',[253 245 157]/255,'Facealpha',0.7)
    hold off
    title(sprintf('Receptor profile (\\gamma = %0.0e)',...
                 gammalist(gammaind(ii))),'FontWeight','normal');
    ylabel('receptor level (%)')
    xlabel('cell surface')
    legend('1 min', '10 min','Location','southeast')
    set(gca,'linewidth',1,'fontsize',12)
    pbaspect([3.5,1,1])
end
nexttile(2,[3 1]);
semilogx(gammalist, speedMAT,'linewidth',1)
xlim([0.0001,0.0316])
xlabel('\gamma')
ylabel('receptor cap speed (\mu m/min)')
set(gca,'linewidth',1,'fontsize',12)
pbaspect([1,1,1])
saveas(gcf,'gamma_timescale.svg')
