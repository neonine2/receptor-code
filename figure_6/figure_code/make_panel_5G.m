function [] = make_panel_5G()

%% Plotting tissue env
load('TissueColormap')
colormap(TissueColormap);
%loading files
load("tissue_300by900_szopt.mat","fconc","fenv")

xlim = 1:0.1:15;
ylim = 474:0.1:488;
coord = combvec(xlim,ylim)';
env = reshape(fconc(coord(:,1),coord(:,2)),length(xlim),length(ylim))';
imagesc(env) 
line([101,130]',[12,12]','Color','white','Linewidth',3) % 3um scalebar
set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal');
colorbar('location','southoutside','fontsize',15,'Direction','reverse')
pbaspect([1,1,1])
saveas(gca,"panel_5G.svg");

end