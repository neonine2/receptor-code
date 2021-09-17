function [] = make_panel_5H(fname)
%% Plotting trajectories
%loading files
load(fname)
load("tissue_300by900_szopt.mat","fconc","fenv")

%trajectory
colormap('gray')
tiledlayout(3,2,'TileSpacing','compact','Padding','none')
vec = [2, 7, 23, 85, 89, 121];
for jj = 1:6
    ii = vec(jj);
    unif = posUnif(1:2:60*2,:,ii);
    opt = posScheme(1:2:60*2,:,ii);
    startpt = opt(1,:);
    
    nexttile;
    xcoord = -5:0.1:10;
    ycoord = -10:0.1:10;
    coord = startpt + combvec(xcoord,ycoord)';
    env = reshape(fconc(coord(:,1),coord(:,2)),length(xcoord),length(ycoord))';
    imagesc(env)

    hold on;
    plot((unif(:,1))*10, (unif(:,2)-startpt(2)+10)*10-3,...
        'Color',[0 0.4470 0.7410]*1.12,...
        'LineWidth',1,...
        'MarkerFaceColor',[0.5,0.5,0.5]);
    plot((opt(:,1))*10, (opt(:,2)-startpt(2)+10)*10-3,...
        'Color',[0.8500 0.3250 0.0980]*1.12,...
        'LineWidth',1,...
        'MarkerFaceColor',[0.5,0.5,0.5]);
    line([111,140]',[15,15]','Color','white') % 3um scalebar
    set(gca,'Xticklabel',[],'Yticklabel',[], 'YDir','normal')
    pbaspect([1,1,1])
    hold off
end
saveas(gca, "panel_5H.svg")
end