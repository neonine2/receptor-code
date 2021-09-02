function [] = make_panel_1A(fnamelist)
%% tissue
load(fnamelist(2),'fconc','fenv')
load('TissueColormap')
colormap(TissueColormap);
[xmin,xmax,ymin,ymax] = deal(1,250,1,250); %um
posmat = combvec(xmin:xmax,ymin:ymax)';
env = reshape(fconc(posmat(:,1),posmat(:,2)),xmax,ymax)';
tiledlayout(1,2)
nexttile
imagesc(env);
set(gca,'YTickLabel',[],'XTickLabel',[],'xtick',[],'ytick',[],'LineWidth',1.5);
pbaspect([1,1,1])
hold on;
line([186,235]',[235,235]','Color','white','Linewidth',3)
MINCOLOR = min(env,[],'all');
MAXCOLOR = max(env,[],'all');
colorbar(gca,'Ticks',[0,1,2],'Ticklabels',["0","1","2"],...
    'location','southoutside','fontsize',14,'Linewidth',0.5); 
caxis([MINCOLOR MAXCOLOR]);
rectangle('Position',[31 31 40 40],'linewidth',2,'EdgeColor',[12,12,12]/17)
hold off;

nexttile
imagesc(imresize(env(31:70,31:70),2.3));
hold on;
line([22,44]'+42,[44,44]'+42,'Color','white','Linewidth',3)
set(gca,'YTickLabel',[],'XTickLabel',[],'xtick',[],'ytick',[],...
    'Linewidth',1.5,'XColor',[12,12,12]/17,'YColor',[12,12,12]/17);
pbaspect([1,1,1])
colorbar(gca,'Ticks',[0,1],'Ticklabels',["0","1"],...
    'location','southoutside','fontsize',14,'Linewidth',0.5);
hold off;
saveas(gcf,'panel_1A_tissue.svg')

%% gradient
load(fnamelist(3),'fconc')
load('GradColormap')
colormap(GradColormap);
[xmin,xmax,ymin,ymax] = deal(1,250,1,250); %um
posmat = combvec(xmin:xmax,ymin:ymax)';
env = reshape(fconc(posmat(:,1),posmat(:,2)),xmax,ymax)';
f = fit((xmin:xmax)',mean(env)','exp1'); % concentration function
env = repmat(f(xmin:xmax)',size(env,1),1);
tiledlayout(1,2)
nexttile
imagesc(env); % shows a piece of the entire environment
set(gca,'YTickLabel',[],'XTickLabel',[],'xtick',[],'ytick',[],'LineWidth',1.5);
pbaspect([1,1,1])
hold on;
line([186,235]',[235,235]','Color','white','Linewidth',3)
MINCOLOR = min(env,[],'all');
MAXCOLOR = max(env,[],'all');
colorbar(gca,'location','southoutside','fontsize',14,'Linewidth',0.5); 
caxis([MINCOLOR MAXCOLOR]);
rectangle('Position',[31 31 40 40],'linewidth',2,'EdgeColor',[12,12,12]/17)
hold off;

nexttile
imagesc(imresize(env(31:70,31:70),2.3));
hold on;
line([22,44]'+42,[44,44]'+42,'Color','white','Linewidth',3)
set(gca,'YTickLabel',[],'XTickLabel',[],'xtick',[],'ytick',[],...
    'Linewidth',1.5,'XColor',[12,12,12]/17,'YColor',[12,12,12]/17);
pbaspect([1,1,1])
colorbar(gca,'location','southoutside','fontsize',14,'Linewidth',0.5);
hold off;
saveas(gcf,'panel_1A_gradient.svg')

%% soil
load(fnamelist(1),'fconc','fenv')
load('SoilColormap')
[xmin,xmax,ymin,ymax] = deal(301,550,301,550); %um
posmat = combvec(xmin:xmax,ymin:ymax)';
envsoil = reshape(fconc(posmat(:,1),posmat(:,2)),xmax-xmin+1,ymax-ymin+1)';

colormap(SoilColormap);
figure(1)
tiledlayout(1,2)
nexttile
imagesc(envsoil); % shows a piece of the entire environment
set(gca,'YTickLabel',[],'XTickLabel',[],'xtick',[],'ytick',[],'LineWidth',1.5);
pbaspect([1,1,1])
hold on;
line([186,235]',[235,235]','Color','white','Linewidth',3)
MINCOLOR = min(envsoil,[],'all');
MAXCOLOR = max(envsoil,[],'all');
colorbar(gca,...
    'location','southoutside','fontsize',14,'Linewidth',0.5); 
caxis([MINCOLOR MAXCOLOR]);
rectangle('Position',[131 141 40 40],'linewidth',2,'EdgeColor',[12,12,12]/17)
hold off;

nexttile
imagesc(imresize(envsoil(131:170,141:180),2.3));
hold on;
line([22,44]'+42,[44,44]'+42,'Color','white','Linewidth',3)
set(gca,'YTickLabel',[],'XTickLabel',[],'xtick',[],'ytick',[],...
    'Linewidth',1.5,'XColor',[12,12,12]/17,'YColor',[12,12,12]/17);
pbaspect([1,1,1])
colorbar(gca,'location','southoutside','fontsize',14,'Linewidth',0.5);
hold off;
saveas(gcf,'panel_1A_soil.svg')

end