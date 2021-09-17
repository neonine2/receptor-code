params = struct('rtot',10,'kd',40,'receptornoise',0.1);
gradientmat = zeros(20,2);
r = linspace(0.01,0.99,20);
for ii = 1:20
    [~,grad1] = totalMI(r(ii),1,compute_exlogx(1,params),params);
    [~,grad2] = totalMI(r(ii),1.05,compute_exlogx(1.05,params),params);
    gradientmat(ii,:) = -[grad1,grad2];
end

params = struct('rtot',40,'kd',40,'receptornoise',0.1);
gradientmat_2 = zeros(20,2);
r = linspace(0.01,0.99,20);
for ii = 1:20
    [~,grad1] = totalMI(r(ii),1,compute_exlogx(1,params),params);
    [~,grad2] = totalMI(r(ii),1.05,compute_exlogx(1.05,params),params);
    gradientmat_2(ii,:) = -[grad1,grad2];
end

params = struct('rtot',200,'kd',40,'receptornoise',0.1);
gradientmat_3 = zeros(20,2);
r = linspace(0.01,0.99,20);
for ii = 1:20
    [~,grad1] = totalMI(r(ii),1,compute_exlogx(1,params),params);
    [~,grad2] = totalMI(r(ii),1.05,compute_exlogx(1.05,params),params);
    gradientmat_3(ii,:) = -[grad1,grad2];
end

%%plotting
tiledlayout(1,3,'tilespacing','compact')
nexttile
plot(r,gradientmat-max(gradientmat,[],'all')+1,'linewidth',2)
yline(0.99935,'-','Optimal','LineWidth',1,'fontsize',16,...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
xline(0.02,'--','LineWidth',1,'color','#0072BD','Alpha',1);
xline(0.98,'--','LineWidth',1,'color','#D95319','Alpha',1);
set(gca,'fontsize',12)
title("$N = 10$","interpreter","latex",'fontsize',18)
ylabel('$\frac{d}{dr_i}I(X;Y)$','interpreter','latex','fontsize',20)
xlabel('$r_i/N$','interpreter','latex','fontsize',20)
pbaspect([1,1,1])
legend("channel 1","channel 2","Location","southwest",...
        'fontsize',13)
legend("boxoff")
ylim([0.99,1])

nexttile
plot(r,gradientmat_2-max(gradientmat_2,[],'all')+1,'linewidth',2)
pbaspect([1,1,1])
yline(0.9938,'-','Optimal','LineWidth',1,'fontsize',16,...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
xline(0.37,'--','LineWidth',1,'color','#0072BD','Alpha',1);
xline(0.63,'--','LineWidth',1,'color','#D95319','Alpha',1);
set(gca,'fontsize',12)
title("$N = 40$","interpreter","latex",'fontsize',18)
xlabel('$r_i/N$','interpreter','latex','fontsize',20)
pbaspect([1,1,1])
ylim([0.98,1])

nexttile
plot(r,gradientmat_3-max(gradientmat_3,[],'all')+1,'linewidth',2)
pbaspect([1,1,1])
yline(0.8975,'-','Optimal','LineWidth',1,'fontsize',16,...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
xline(0.475,'--','LineWidth',1,'color','#0072BD','Alpha',1);
xline(0.522,'--','LineWidth',1,'color','#D95319','Alpha',1);
set(gca,'fontsize',12)
title("$N = 200$","interpreter","latex",'fontsize',18)
xlabel('$r_i/N$','interpreter','latex','fontsize',20)
pbaspect([1,1,1])
ylim([0.85,1])

saveas(gca,"derivative_rtot.svg")