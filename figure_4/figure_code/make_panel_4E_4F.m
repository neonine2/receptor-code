function [] = make_panel_4E_4F(fname)
load(fname)
recC = results.env;
recF = results.f;
recC = recC((160*60+1):60:20400,:);
recF = recF((160*60+1):60:20400,:);
[~,recCmax] = max(recC,[],2);
recCmax = squeeze(recCmax);
[~,recFmax] = max(recF,[],2);
recFmax = squeeze(recFmax);
m = scheme_param.N;

% plotting kymograph and max positions
figure(1)
tiledlayout(2,1,'tilespacing','compact')
c = viridis;
colormap(c);
nexttile
imagesc(flipud(recF'))
colorbar('Ticks',[10,50],'TickLabels',["10","50"])
yticks([1,50,100])
yticklabels({'\pi','0','-\pi'})
set(gca,'fontsize',17)
xticks(30*(1:8))
xticklabels({'30','60','90','120','150','180','210'})
set(gca,'xaxisLocation','top')
pbaspect([4,1,1])

nexttile
hold on
plot(recCmax,'color',[12 12 12]/20,'linewidth',2);
plot(recFmax,'color',[0.9290 0.6940 0.1250],'linewidth',2);
xlim([1,length(recCmax)])
ylim([1,100])
yticks([1,50,100])
yticklabels({'-\pi','0','\pi'})
set(gca,'xticklabels',[])
hold off
set(gca,'fontsize',17)
pbaspect([4,1,1])
box on

saveas(gca,"panel_4D.svg")

%% panel F snapshots
bdd = 18;
xlimit = -bdd:0.1:bdd;
figure(2)
tiledlayout(1,4)
load('tissue_300by900_szopt','fconc','fenv');
snapshotindex = [11221,13261,15541,19841];
s = scheme_param.mean_cell_radius*10;
for ii = 1:4
    colormap(gray)
    nexttile
    indextissue = snapshotindex(ii);
    coord = combvec(xlimit,xlimit)' + results.cellp(indextissue,:);
    envtissue = reshape(arrayfun(fconc,coord(:,1),coord(:,2)),...
        length(xlimit),length(xlimit))';
    imagesc(envtissue)
    pbaspect([1,1,1])
    r = results.f(indextissue,:)';
    angle = linspace(0,2*pi*(1-1/m),m);
    cellsurf = 181 + (s+r*1.4).*[cos(angle)',sin(angle)'];
    hold on
    fill(cellsurf(:,1),cellsurf(:,2),[17 17 17]/18,'Linewidth',1)
    rectangle('Position',[181-s,181-s,2*s,2*s],'Curvature',[1,1], ...
        'FaceColor', [17 17 17]/32, 'EdgeColor', 'k','Linewidth',2);
    set(gca,'xticklabels',[],'yticklabels',[],'xtick',[],'ytick',[]);
    hold off
end

saveas(gca,"panel_4E.svg")

%% panel G
figure(3)
histogram(recCmax,30,'FaceColor',[12 12 12]/20,'EdgeColor','none',...
    'Normalization','probability')
hold on
histogram(recFmax,11,'FaceColor',[0.9290 0.6940 0.1250],...
    'EdgeColor','none','Normalization','probability')
hold off
ylim([0,0.3])
xlabel("Membrane position (rad)")
ylabel("Probability")
legend("Ligand peak","Receptor peak",'fontsize',10)
legend('boxoff')
pbaspect([1,2,1])
box off
set(gca,'fontsize',12,'xticklabels',["-\pi","0","\pi"])
saveas(gca,"panel_4G.svg")

end
