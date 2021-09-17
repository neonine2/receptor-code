function [] = make_panel_5B()

%% Plotting tissue env
load('TissueColormap')
colormap(TissueColormap);
%loading files
load("tissue_300by900_szopt.mat","fconc","fenv")

xlim = 1:0.1:70;
ylim = 501:0.1:570;
coord = combvec(xlim,ylim)';
env = reshape(fconc(coord(:,1),coord(:,2)),length(xlim),length(ylim))';
imagesc(env) 
line([501,600]',[40,40]','Color','white','Linewidth',3) %10um scalebar
set(gca,'YTickLabel',[],'XTickLabel',[],'LineWidth',2,'YDir','normal');
colorbar('location','southoutside','fontsize',15,'Direction','reverse')
pbaspect([1,1,1])

saveas(gca,"panel_5B.svg")
end