function [] = make_panel_1D(fnamelist)

%% loading relative efficacy
%tissue and grad
load('colorpalette')
figurec(1)
load(fnamelist(1))
load(fnamelist(2))
nkecm = 9;
kecmlist = logspace(-4,-2,nkecm)*12/520;

tiledlayout(1,2)
%soil
load(fnamelist(3))
nexttile
varlist = [logspace(-3,1.5,10), 80];
semilogx(varlist(5:end),rel_eff_soil(5:end),'linewidth',1,'color',hex2rgb(hexcolor{1}))
ylim([40,220])
xlim([varlist(5),varlist(end)])
set(gca,'fontsize',16)
pbaspect([2.8,1,1])

nexttile
semilogx(kecmlist,rel_eff_tissue,'color',hex2rgb(hexcolor{2}),...
    'linewidth',1)
hold on
semilogx(kecmlist,rel_eff_grad,'color',hex2rgb(hexcolor{3}),...
            'linewidth',1)
hold off 
xlim([min(kecmlist),max(kecmlist)])
% xtix = get(gca, 'XTick');
% set(gca, 'XTick',xtix, 'XTickLabel',xtix*12/520);
set(gca,'fontsize',16)
pbaspect([2.8,1,1])

end