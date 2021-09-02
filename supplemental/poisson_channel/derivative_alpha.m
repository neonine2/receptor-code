% 2 channel run
clear
m = 2;
epsilon = 0.01./m;
Aeq = ones(1,m);
beq = 1;
lb = zeros(1,m);
ub = ones(1,m);
params = struct('rtot',20,'kd',40,'receptornoise',0.1);
options = optimoptions('fmincon','Algorithm','sqp',...
                       'SpecifyObjectiveGradient',true,...
                       'Display','off'); 
                   
% make plot of change vs. mu
nalpha = 100;
mu = [0.1,1.05*0.1];
alphavec = logspace(-3,0,nalpha);
nkd = 100;
kdvec = logspace(0,3,nkd);
paramlist = combvec(alphavec,kdvec)';
optR = zeros(nalpha*nkd,1);
approxMIslope = zeros(nalpha*nkd,1);
parfor ii = 1:length(paramlist)
    newparams = params;
    newparams.receptornoise = paramlist(ii,1);
    newparams.kd = paramlist(ii,2);
    exlogx = compute_exlogx(mu,newparams);
    objective = @(x) totalMI(x,mu,exlogx,newparams);
    x0 = mu'./sum(mu);
    [x,fval] = fmincon(objective,x0,[],[],Aeq,beq,lb,ub,[],options);
    optR(ii) = x(1);
    approxMIslope(ii) = approx_grad(mu(1),newparams);
end


c = viridis;
colormap(c);
imagesc(kdvec, alphavec, reshape(hilldiff,nalpha,nkd));
set(gca, 'fontsize',12,'XScale','log','YScale','log')
xlabel('$K_d$','interpreter','latex','fontsize',24)
ylabel('$\alpha$','interpreter','latex','fontsize',24)
ylim(10.^[-3,0])
xlim(10.^[0,3])
colorbar()
pbaspect([1,1,1])

c = viridis;
colormap(c);
imagesc(kdvec, alphavec, reshape(optR,nalpha,nkd));
set(gca, 'fontsize',12,'XScale','log','YScale','log')
xlabel('$K_d$','interpreter','latex','fontsize',24)
ylabel('$\alpha$','interpreter','latex','fontsize',24)
ylim(10.^[-3,0])
xlim(10.^[0,3])
colorbar()
caxis([0,0.5])
pbaspect([1,1,1])

