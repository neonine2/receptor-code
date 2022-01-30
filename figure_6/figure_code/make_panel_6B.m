function [] = make_panel_6B(filename,gradfilename)

%% Plotting tissue env
load('TissueColormap')
colormap(TissueColormap);
%loading files
load(filename,"fconc","fenv")

xlim = 1:0.1:70;
ylimit = 501:0.1:570;
coord = combvec(xlim,ylimit)';
env = reshape(fconc(coord(:,1),coord(:,2)),length(xlim),length(ylimit))';

tiledlayout(2,1,'TileSpacing','Compact');
nexttile;
imagesc(env) 
line([501,600]',[40,40]','Color','white','Linewidth',3) %10um scalebar
set(gca,'YTickLabel',[],'XTickLabel',[],'LineWidth',2,'YDir','normal');
colorbar('location','southoutside','fontsize',15,'Direction','reverse')
pbaspect([1,1,1])

nexttile;
meanprofile = mean(env);
f = fit(xlim',meanprofile','exp1');
plot(xlim,meanprofile,'r',xlim,f(xlim),'k')
ylim([min(meanprofile),max(meanprofile)])
box off
pbaspect([1,0.2,1])

saveas(gca,"panel_6B.svg")
end