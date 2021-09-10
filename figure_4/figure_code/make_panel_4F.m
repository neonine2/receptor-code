function [] = make_panel_4G(fname)

load('colorpalette')
%% loading hlist data
load(fname(1),'relative_eff','recstat');
meankfb_tissue = squeeze(mean(recstat(:,1,:)));
rel_eff_tissue = relative_eff;
load(fname(2),'relative_eff','recstat');
meankfb_soil = squeeze(mean(recstat(:,1,:)));
rel_eff_soil = relative_eff;

% plotting h para
tiledlayout(2,1,'tilespacing','compact')
nexttile
semilogx(meankfb_soil,rel_eff_soil,'Color',hex2rgb(hexcolor{1}),'linewidth',2)
hold on
semilogx(meankfb_tissue,rel_eff_tissue,'Color',hex2rgb(hexcolor{2}),'linewidth',2)
hold off
ylim([0.1,0.8])
xlim([10^(-4),10^(-2.5)])
pbaspect([2,1,1])
set(gca,'fontsize',16,'linewidth',1.5)

%% loading kofflist data
load(fname(3),'relative_eff');
rel_eff_tissue = relative_eff;
load(fname(4),'relative_eff','kofflist');
rel_eff_soil = relative_eff;

% plotting koff para
nexttile
semilogx(kofflist,rel_eff_soil,'Color',hex2rgb(hexcolor{1}),'linewidth',2)
hold on
semilogx(kofflist,rel_eff_tissue,'Color',hex2rgb(hexcolor{2}),'linewidth',2)
hold off
xlim([1e-2,10^(-0.6)])
ylim([0,0.8])
pbaspect([2,1,1])
set(gca,'fontsize',16,'linewidth',1.5)

saveas(gca,"panel_4G.svg")
end