function [] = make_panel_1B(fnamelist)
%% loading data files
load('colorpalette')
cellmodel = 2; % 10um cell

nfile = length(fnamelist);
%% panel B-i)
figure(1)
tiledlayout(1,3,'TileSpacing','none')

%simple gradient
load(fnamelist(3),'envmean','optr')
nexttile
yyaxis left
plot(optr(11,:,cellmodel)*100,'Linewidth',2,'color',hexcolor{3})
% sum(maxk(optr(11,:,cellmodel),12))
% -sum(optr(11,:,cellmodel).*log(optr(1+10,:,cellmodel)),'omitnan')
ylim([0,11])
yyaxis right
plot(envmean(1+10,:,cellmodel),'Linewidth',2,'color',[1,1,1]/2)
pbaspect([1,1,1])
ylim([0,1.2])
set(gca,'YColor','none')
box off
set(gca,'Linewidth',1.5,'fontsize',13,'YTick',[0,0.5,1],...
    'XTickLabel',{'-\pi','0','\pi'});

%soil
load(fnamelist(1),'envmean','optr')
nexttile
yyaxis left
plot(optr(67,:,cellmodel)*100,'Linewidth',2,'color',hexcolor{1})
ylim([0,11])
% sum(maxk(optr(67,:,cellmodel),12))
% -sum(optr(67,:,cellmodel).*log(optr(67,:,cellmodel)),'omitnan')
set(gca,'YColor','none')
yyaxis right
plot(envmean(67,:,cellmodel),'Linewidth',2,'color',[1,1,1]/2)
pbaspect([1,1,1])
ylim([0,1.2])
box off
set(gca,'YColor','none')
set(gca,'Linewidth',1.5,'fontsize',13,'YTick',[0,0.5,1],...
    'XTickLabel',{'-\pi','0','\pi'});

%tissue
load(fnamelist(2),'envmean','optr')
nexttile
yyaxis left
plot(optr(81,:,cellmodel)*100,'Linewidth',2,'color',hexcolor{2})
ylim([0,11])
% sum(maxk(optr(81,:,cellmodel),12))
% -sum(optr(81,:,cellmodel).*log(optr(81,:,cellmodel)),'omitnan')
set(gca,'YColor','none')
yyaxis right
plot(envmean(81,:,cellmodel),'Linewidth',2,'color',[1,1,1]/2)
pbaspect([1,1,1])
ylim([0,1.2])
box off
set(gca,'Linewidth',1.5,'fontsize',13,'YTick',[0,0.5,1],...
    'XTickLabel',{'-\pi','0','\pi'});
saveas(gca,"panel_1B_lineplot.svg")

%% panel B-ii)
% for each environment, compute entropy of receptor placement
proportioncell = cell(1,length(fnamelist));
for ii = 1:nfile
    load(fnamelist(ii),'envmean','optr')
    [nloc, m] = size(envmean,[1,2]);
    optr = permute(optr(:,:,cellmodel),[1 3 2]);
    optr = reshape(optr,[],m,1);
    proportion = zeros(nloc,1);
    for jj = 1:nloc
        rvec = optr(jj,:)./sum(optr(jj,:));
        rvec = rvec(rvec~=0);
        proportion(jj) = -sum(rvec.*log(rvec));
    end
    proportioncell{ii} = proportion;
end

% plotting histogram
figure(2)
tiledlayout(1,2)
nexttile
hold on
for ii = 1:nfile
    % put soil in front of tissue for better visual
    if ii == 1
        ii = 2;
    elseif ii == 2
        ii = 1;
    end
    histogram(proportioncell{ii},'BinWidth',0.1,...
        'Normalization','probability',...
            'FaceAlpha',0.8,'Facecolor',hexcolor{ii})
end
pbaspect([1.5,1,1])
ylim([0,0.12])
ytix = get(gca, 'YTick');
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100);
set(gca,'Linewidth',1,'fontsize',16);
xlim([0,-log(1/m)])
box off
hold off

%% panel B-iii)
nexttile
IcMAT = zeros(300,3);
IrMAT = zeros(300,3);

for ii = 1:3
    load(fnamelist(ii),'envmean','optr')
    [~, Ic] = max(envmean(1:300,:,cellmodel),[],2);
    [~, Ir] = max(optr(1:300,:,cellmodel),[],2);
    IcMAT(:,ii) = Ic./m.*2.*pi-pi;
    IrMAT(:,ii) = Ir./m.*2.*pi-pi;
end

% plotting scatterplot
hold on
for ii = 1:nfile
    scatter(IcMAT(:,ii),IrMAT(:,ii),'filled','MarkerFaceAlpha',0.2,...
        'MarkerFaceColor', hexcolor{ii})
end
hold off
set(gca, 'XTick',[-pi,0,pi], 'XTickLabel',{'-\pi','0','\pi'},'Xlim',[-pi,pi]);
set(gca, 'YTick',[-pi,0,pi], 'YTickLabel',{'-\pi','0','\pi'},'Ylim',[-pi,pi]);
set(gca,'Linewidth',1,'fontsize',16)
pbaspect([1,1,1])
saveas(gca,"panel_1B.svg")

end

