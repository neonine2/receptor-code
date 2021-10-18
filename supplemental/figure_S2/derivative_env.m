% 2 channel run
clear all;
close all;

m = 2;
epsilon = 0.01./m;
Aeq = ones(1,m);
beq = 1;
lb = zeros(1,m);
ub = ones(1,m);
params = struct('rtot',40,'kd',40,'receptornoise',0.1);
options = optimoptions('fmincon','Algorithm','sqp',...
                       'SpecifyObjectiveGradient',true,...
                       'Display','off'); 
                   
% make plot of change vs. mu
nmean = 100;
muvec = logspace(-1,1.7,nmean);
optR = zeros(nmean,1);
optMI = zeros(nmean,1);
exvec = zeros(nmean,1);
approxMIslope = zeros(nmean,1);
for ii = 1:nmean
    mu = [muvec(ii),1.05*muvec(ii)];
    exlogx = compute_exlogx(mu,params);
    objective = @(x) totalMI(x,mu,exlogx,params);
    x0 = mu;
    x0 = x0/sum(x0);
    [x,fval] = fmincon(objective,x0,[],[],Aeq,beq,lb,ub,[],options);
%     for jj = 1:500
%         mivec(jj) = totalMI([jj/500,1-jj/500], mu, exlogx, params);
%     end
%     [fval, x] = min(mivec);
    optR(ii) = x(1);
    optMI(ii) = -fval;
    [grad,ex] = approx_grad(muvec(ii),params);
    approxMIslope(ii) = grad;
    exvec(ii) = ex/params.rtot;
end

tiledlayout(1,2)
left_color = [0 0 0];
right_color = [0 0.4470 0.7410];
set(figure,'defaultAxesColorOrder',[left_color; right_color]);
nexttile
yyaxis left
plot(exvec,optR,'k','linewidth',1);
yline(0.5,'--r','linewidth',1)
set(gca,'fontsize',16, 'XScale', 'log')
ylabel('$r^*_1/N$','interpreter','latex','fontsize',23)
yyaxis right
plot(exvec,approxMIslope,'linewidth',1)
[~,ind] = max(approxMIslope);
xline(exvec(ind),'--', 'linewidth',1)
ylabel('$\partial I_0/\partial r_1$','interpreter','latex','fontsize',23,...
    'Rotation',-90,'VerticalAlignment','bottom')
xlabel('$E[X_1]$','interpreter','latex','fontsize',23)
% xlim([10^-0.5,10^1.5])
patch([exvec(ind) exvec(ind) max(xlim) max(xlim)],...
    [0 0.2 0.2 0],...
    [0.9290 0.6940 0.1250],'FaceAlpha',.3,'EdgeColor','none')
text(0.2,0.012,"$1.05 \times E[X_1] = E[X_2]$",'interpreter','latex',...
    'fontsize',18)
pbaspect([1,1,1])

nexttile
pvec = linspace(0,1,100);
plot(pvec,-pvec.*log(pvec)-(1-pvec).*log(1-pvec),'linewidth',1)
xlabel('Success probability p','fontsize',23)
ylabel("H(p)",'fontsize',23,'interpreter','latex')
pbaspect([0.5,1,1])
set(gca,'fontsize',16)
saveas(gca,"derivative_env.svg")
