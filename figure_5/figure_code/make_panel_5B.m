function [envtissue,envsoil] = make_panel_5B(fname)
%% Panel B
cellmodel = 2;
bdd = 22;
xlimit = -bdd:0.1:bdd;


indextissue = 19;
load('tissue_300by900_szopt','centerlist','fconc')
load(fname(1))
schemerT = schemer;
posindT = posind;
coord = combvec(xlimit,xlimit)' + centerlist(posindT(indextissue),:,cellmodel);
envtissue = reshape(arrayfun(fconc,coord(:,1),coord(:,2)),...
    length(xlimit),length(xlimit))';


indexsoil = 20;
load('soil_var_2_szopt','centerlist','fconc')
load(fname(2))
schemerS = schemer;
posindS = posind;
coord = combvec(xlimit,xlimit)' + centerlist(posindS(indexsoil),:,cellmodel);
envsoil = reshape(arrayfun(fconc,coord(:,1),coord(:,2)),...
    length(xlimit),length(xlimit))';


%% plotting
colormap(gray)
tiledlayout(2,2,'TileSpacing','compact')
nexttile
imagesc(envtissue)
pbaspect([1,1,1])
m = scheme_param.N;
r = ones(m,1)./m;
angle = linspace(0,2*pi*(1-1/m),m);
s = scheme_param.mean_cell_radius*10;
cellsurf = 221 + (s+r*m*10).*[cos(angle)',sin(angle)'];
hold on
fill(cellsurf(:,1),cellsurf(:,2),[17 17 17]/18,'Linewidth',1)
rectangle('Position',[221-s,221-s,2*s,2*s],'Curvature',[1,1], ...
    'FaceColor', [17 17 17]/32, 'EdgeColor', 'k','Linewidth',2);
set(gca,'xticklabels',[],'yticklabels',[],'xtick',[],'ytick',[],...
    'linewidth',2);


nexttile
imagesc(envtissue)
pbaspect([1,1,1])
m = scheme_param.N;
r = schemerT(indextissue,:,24)';
angle = linspace(0,2*pi*(1-1/m),m);
s = scheme_param.mean_cell_radius*10;
cellsurf = 221 + (s+r*m*5).*[cos(angle)',sin(angle)'];
hold on
fill(cellsurf(:,1),cellsurf(:,2),[17 17 17]/18,'Linewidth',1)
rectangle('Position',[221-s,221-s,2*s,2*s],'Curvature',[1,1], ...
    'FaceColor', [17 17 17]/32, 'EdgeColor', 'k','Linewidth',2);
set(gca,'xticklabels',[],'yticklabels',[],'xtick',[],'ytick',[],...
    'linewidth',2);
colorbar('Ticks',[0.5,1,1.5],'TickLabels',["0.5","1","1.5"],...
    'fontsize',12)

nexttile
imagesc(envsoil)
pbaspect([1,1,1])
m = scheme_param.N;
r = ones(m,1)./m;
angle = linspace(0,2*pi*(1-1/m),m);
s = scheme_param.mean_cell_radius*10;
cellsurf = 221 + (s+r*m*10).*[cos(angle)',sin(angle)'];
hold on
fill(cellsurf(:,1),cellsurf(:,2),[17 17 17]/18,'Linewidth',1)
rectangle('Position',[221-s,221-s,2*s,2*s],'Curvature',[1,1], ...
    'FaceColor', [17 17 17]/32, 'EdgeColor', 'k','Linewidth',2);
set(gca,'xticklabels',[],'yticklabels',[],'xtick',[],'ytick',[],...
    'linewidth',2);

nexttile
imagesc(envsoil)
pbaspect([1,1,1])
m = scheme_param.N;
r = schemerS(indexsoil,:,24)';
angle = linspace(0,2*pi*(1-1/m),m);
s = scheme_param.mean_cell_radius*10;
cellsurf = 221 + (s+r*m*5).*[cos(angle)',sin(angle)'];
hold on
fill(cellsurf(:,1),cellsurf(:,2),[17 17 17]/18,'Linewidth',1)
rectangle('Position',[221-s,221-s,2*s,2*s],'Curvature',[1,1], ...
    'FaceColor', [17 17 17]/32, 'EdgeColor', 'k','Linewidth',2);
set(gca,'xticklabels',[],'yticklabels',[],'xtick',[],'ytick',[],...
    'linewidth',2);
colorbar('Ticks',[1,2,3],'TickLabels',["1","2","3"],'fontsize',12)

saveas(gca,"panel_5B.svg")
end