function [] = make_panel_6C(fname)
%% Plotting trajectories
%loading files
load(fname)
load("tissue_300by900_szopt.mat","fconc","fenv")

%trajectory
colormap('gray')
tiledlayout(3,2,'TileSpacing','compact','Padding','none')
vec = [2, 22, 48, 78, 135, 146];
for jj = 1:6
    ii = vec(jj);
    unif = posUnif(1:60*2,:,ii);
    opt = posScheme(1:60*2,:,ii);
    startpt = opt(1,:);

    nexttile;
    xcoord = -45:0.1:10;
    ycoord = -35:0.1:35;
    coord = startpt + combvec(xcoord,ycoord)';
    env = reshape(fconc(coord(:,1),coord(:,2)),length(xcoord),length(ycoord))';
    imagesc(env)

    hold on;
    lasti = find(unif(:,1)==0,1,'first')-1; %finding final coordinate
    if isempty(lasti)
        lasti = 60*2;
    end
    plot((unif(1:lasti,1))*10, (unif(1:lasti,2)-startpt(2)+36)*10,...
        'Color',[0 0.4470 0.7410]*1.12,...
        'LineWidth',2,...
        'MarkerFaceColor',[0.5,0.5,0.5]);
    lasti = find(opt(:,1)==0,1,'first')-1;
    if isempty(lasti)
        lasti = 60*2;
    end
    plot((opt(1:lasti,1))*10, (opt(1:lasti,2)-startpt(2)+36)*10,...
        'Color',[0.8500 0.3250 0.0980]*1.12,...
        'LineWidth',2,...
        'MarkerFaceColor',[0.5,0.5,0.5]);
    line([401,500]',[50,50]','Color','white','Linewidth',2) %10um scalebar
    set(gca,'Xticklabel',[],'Yticklabel',[], 'YDir','normal',...
        'linewidth',1)
    pbaspect([1,1,1])
    hold off
end
brighten(0);

saveas(gca,"panel_6C.svg")
end