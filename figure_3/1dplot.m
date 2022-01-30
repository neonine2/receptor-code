load(which("tissue_300by900_szopt_perturb_10um.mat"))

%compute proportion of receptors without a small neighbourhood of peak
m = size(perturbr,2);
nloc = size(perturbr,1);
nhbd = 0;

for ii = 1:nparam
    [~,I] = max(perturbr(:,:,ii),[],2);
    for jj = 1:nloc
        shiftr = circshift(perturbr(jj,:,ii),m/2-I(jj));
        peakAmount(jj,ii) = sum(shiftr((m/2-nhbd):(m/2+nhbd)));
    end
end
nflatten = length(flattenparam);
nshift = length(shiftparam);
peakAmount = reshape(mean(peakAmount), nshift, nflatten);
perturbMI = reshape(mean(perturbMI), nshift, nflatten);
peakAmounttissue = peakAmount(shiftparam == 0,:);
perturbMI = perturbMI(shiftparam == 0,:);

load(which("tissue_300by900_szopt.mat"),'unifMI')
unifMI = mean(unifMI(:,2));
eff = (perturbMI - unifMI);
efftissue = eff./eff(1)


load(which("soil_var_2_szopt_perturb_10um.mat"))

%compute proportion of receptors without a small neighbourhood of peak
m = size(perturbr,2);
nloc = size(perturbr,1);
nhbd = 0;

for ii = 1:nparam
    [~,I] = max(perturbr(:,:,ii),[],2);
    for jj = 1:nloc
        shiftr = circshift(perturbr(jj,:,ii),m/2-I(jj));
        peakAmount(jj,ii) = sum(shiftr((m/2-nhbd):(m/2+nhbd)));
    end
end
nflatten = length(flattenparam);
nshift = length(shiftparam);
peakAmount = reshape(mean(peakAmount), nshift, nflatten);
perturbMI = reshape(mean(perturbMI), nshift, nflatten);
peakAmountsoil = peakAmount(shiftparam == 0,:);
perturbMI = perturbMI(shiftparam == 0,:);

load(which("soil_var_2_szopt.mat"),'unifMI')
unifMI = mean(unifMI(:,2));
eff = (perturbMI - unifMI);
effsoil = eff./eff(1)

tiledlayout(1,2)
n=50;
mu = logspace(-3,1,n);
ent = zeros(n,1);
for ii = 1:n
    ent(ii) = compute_entropy(mu(ii));
end
nexttile
plot(mu,ent,'linewidth',1)
set(gca,'fontsize',14,'linewidth',1);
pbaspect([1.5,1,1])
xlabel('\lambda')
ylabel('entropy')

nexttile
plot(peakAmountsoil/peakAmountsoil(1), effsoil,'linewidth',1);
hold on
plot(peakAmounttissue/peakAmounttissue(1), efftissue,'linewidth',1);
hold off
set(gca,'fontsize',14,'linewidth',1);
pbaspect([1.5,1,1])
legend('Tissue','Soil','Location','southeast','box','off')
xlabel('proportion of original receptor count')
ylabel('proportion of original efficacy')


