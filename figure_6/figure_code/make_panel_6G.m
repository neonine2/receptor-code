function [] = make_panel_6G(filename)

%% Plotting tissue env
load('TissueColormap')
colormap(TissueColormap);
%loading files
load(filename,"fconc","fenv")

xlimit = 1:0.1:15;
ylimit = 474:0.1:488;
coord = combvec(xlimit,ylimit)';
env = reshape(fconc(coord(:,1),coord(:,2)),length(xlimit),length(ylimit))';

tiledlayout(2,1,'TileSpacing','Compact');
nexttile;
imagesc(env) 
line([101,130]',[12,12]','Color','white','Linewidth',3) % 3um scalebar
set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal');
colorbar('location','southoutside','fontsize',15,'Direction','reverse')
pbaspect([1,1,1])

nexttile;
meanprofile = mean(env);
f = fit(xlimit',meanprofile','exp1');
plot(xlimit,meanprofile,'r',xlimit,f(xlimit),'k')
ylim([min(meanprofile),max(meanprofile)])
xlim([1,15])
box off
pbaspect([1,0.2,1])

saveas(gca,"panel_6G.svg");
end