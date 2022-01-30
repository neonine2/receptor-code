%% generating data
clear all;
close all;
parameter_scan("localization", "parameter_scan_localization_feedback");
parameter_scan("retention", "parameter_scan_retention_feedback");

%% organizing data

% scattered interpolation for to create phase diagram over parameter space
clear
load('parameter_scan_localization_feedback','scheme_succ',"scheme_stat","nval")
Finterp = scatteredInterpolant(scheme_stat(:,2),scheme_stat(:,1),scheme_succ);
kofflist = logspace(-2,-0.7,nval*8);
kfblist = logspace(-4,-2,nval*10);
paramlist = combvec(kofflist,kfblist)';
successMAT = reshape(Finterp(paramlist(:,1),paramlist(:,2)),nval*8,nval*10);

% scattered interpolation for to create phase diagram over parameter space
load('parameter_scan_retention_feedback','scheme_stat',"scheme_err","nval")
Finterp = scatteredInterpolant(scheme_stat(:,2),scheme_stat(:,1),scheme_err);
errorMAT = reshape(Finterp(paramlist(:,1),paramlist(:,2)),nval*8,nval*10);

%% plotting data
c = viridis;
colormap(c);
imagesc(kfblist,kofflist,successMAT)
colorbar('location','northoutside')
set(gca, 'fontsize',16,'XScale','log','YScale','log')
xlim([min(kfblist),max(kfblist)])
ylim([min(kofflist),max(kofflist)])
pbaspect([6,2.7,1])
saveas(gca, "parameter_scan_localization.svg")

c = viridis;
colormap(c);
imagesc(kfblist,kofflist,log10(errorMAT))
cbh = colorbar('location','northoutside');
cbh.Ticks = log10([2.5,5,10,20,40]);
cbh.TickLabels = 10.^(cbh.Ticks);
xlim([min(kfblist),max(kfblist)])
ylim([min(kofflist),max(kofflist)])
set(gca, 'fontsize',16,'XScale','log','YScale','log')
pbaspect([6,2.7,1])
saveas(gca, "parameter_scan_retention.svg")
