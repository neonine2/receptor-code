function [] = make_panel_6G(filename)

%% Plotting tissue env
load('TissueColormap')
colormap(TissueColormap);
%loading files
load(filename,"fconc","fenv")

xlim = 1:0.1:15;
ylim = 474:0.1:488;
coord = combvec(xlim,ylim)';
env = reshape(fconc(coord(:,1),coord(:,2)),length(xlim),length(ylim))';
imagesc(env) 
line([101,130]',[12,12]','Color','white','Linewidth',3) % 3um scalebar
set(gca,'YTickLabel',[],'XTickLabel',[],'YDir','normal');
colorbar('location','southoutside','fontsize',15,'Direction','reverse')
pbaspect([1,1,1])
saveas(gca,"panel_6G.svg");

end