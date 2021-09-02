% 2 channel run
clear
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
nmean = 20;
muvec = linspace(1,1.1,nmean);
nfac = 5;
facvec = logspace(-0.6,1.6,nfac);
optR = zeros(nmean,1);
for jj = 1:nfac
    for ii = 1:nmean
        mu = [1,muvec(ii)]*facvec(jj);
        exlogx = compute_exlogx(mu,params);
        objective = @(x) totalMI(x,mu,exlogx,params);
        x0 = mu'./sum(mu);
        [x,fval] = fmincon(objective,x0,[],[],Aeq,beq,lb,ub,[],options);
        optR(ii,jj) = x(1);
    end
end

hold on
for k = 1:nfac
    plot((muvec-1)./1,optR(:,k),'linewidth',1.2,...
        'Color',[0.8 0.8 0.8]-0.15*(k-1));
end
hold off
set(gca,'fontsize',16)
lgd = legend(num2str(round(facvec(1),1)),...
        num2str(round(facvec(2),1)),...
        num2str(round(facvec(3),1)),...
        num2str(round(facvec(4),1)),...
        num2str(round(facvec(5),1)),'fontsize',15);
title(lgd,'$E[C_1]$', 'interpreter', 'latex','fontsize',17)
legend('boxoff')
ylabel('$r^*_1/N$','interpreter','latex','fontsize',23)
xlabel('$(E[C_2]-E[C_1])/E[C_1]$','interpreter','latex','fontsize',23)
yline(0.5,'--r','uniform receptor','LabelVerticalAlignment','middle',...
    'linewidth',0.8,'fontsize',14)
xlim([0,0.1])
pbaspect([1,1,1])
box on
saveas(gca,"derivative_drange.svg")

