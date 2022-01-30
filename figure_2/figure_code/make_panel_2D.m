function [] = make_panel_2D(tissue_filelist,soil_filelist,mode)

%% loading relative efficacy
%tissue and grad
load('colorpalette')
figurec(1)
npar = length(tissue_filelist);
eff_grad = zeros(npar,1);
for ii = 1:npar
    load(strcat(tissue_filelist(ii),"_grad_opt.mat"),"optMI","unifMI")
    eff_grad(ii) = compute_efficacy(mean(optMI), mean(unifMI),"mode",mode);
end
eff_tissue = zeros(npar,1);
for ii = 1:npar
    load(strcat(tissue_filelist(ii),"_opt.mat"),"optMI","unifMI")
    eff_tissue(ii) = compute_efficacy(mean(optMI), mean(unifMI),"mode",mode);
end
nkecm = 9;
kecmlist = logspace(-4,-2,nkecm)*12/520;

tiledlayout(1,2)
%soil
npar = length(soil_filelist);
eff_soil = zeros(npar,1);
for ii = 1:npar
    load(strcat(soil_filelist(ii),"_opt.mat"),"optMI","unifMI")
    eff_soil(ii) = compute_efficacy(mean(optMI), mean(unifMI),"mode",mode);
end

nexttile
varlist = [logspace(-3,1.5,10), 80];
semilogx(varlist(5:end),eff_soil(5:end),...
    'linewidth',1,'color',hex2rgb(hexcolor{1}))
% ylim([0,0.7])
xlim([varlist(5),varlist(end)])
set(gca,'Linewidth',1,'fontsize',16)
pbaspect([2.8,1,1])

nexttile
semilogx(kecmlist,eff_tissue,'color',hex2rgb(hexcolor{2}),...
    'linewidth',1)
hold on
semilogx(kecmlist,eff_grad,'color',hex2rgb(hexcolor{3}),...
            'linewidth',1)
hold off 
xlim([min(kecmlist),max(kecmlist)])
% ylim([0.2,0.5])
% xtix = get(gca, 'XTick');
% set(gca, 'XTick',xtix, 'XTickLabel',xtix*12/520);
set(gca,'Linewidth',1,'fontsize',16)
pbaspect([2.8,1,1])
saveas(gca,"panel_1D.png")

end