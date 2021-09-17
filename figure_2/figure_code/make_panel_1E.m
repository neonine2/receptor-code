function [] = make_panel_1E(fnamelist)
load('colorpalette','hexcolor')
nenv = length(fnamelist);
load(fnamelist(2),'envmean')
[nloc,m,ncell] = size(envmean);
sparsemat = zeros(nloc*ncell,nenv);
avgmat = zeros(nloc*ncell,nenv);
deltamat = zeros(nloc*ncell,nenv);
% maxsparsity = -log(1/m);

for ii = 1:nenv
    load(fnamelist(ii),'envmean','optMI','unifMI')
    deltaI = (optMI-unifMI)./unifMI .* 100;
    %remove -ve values
    deltaI(deltaI<0) = 0/0; % make negative nan
    sparsity = 1-squeeze((vecnorm(envmean,1,2)./m)./sqrt((vecnorm(envmean,2,2).^2)/m));
%     sparsity = 1-sum(-envmean./sum(envmean,2).*log(envmean./sum(envmean,2)),2)/maxsparsity;
    avgconc = squeeze(mean(envmean,2));
    deltamat(:,ii) = deltaI(:);
    avgmat(:,ii) = avgconc(:);
    sparsemat(:,ii) = sparsity(:);
end

% making scatterplot
order = [1,3,2];
for jj = 1:nenv
ii = order(jj);
jitterAmount = 3;
jitterValuesX = rand(size(avgmat(:,ii)))*jitterAmount;   % +/-jitterAmount max
scatter3(avgmat(:,ii).*jitterValuesX,...
         sparsemat(:,ii),...
         deltamat(:,ii),45,hex2rgb(hexcolor{ii}),...
         'filled','MarkerFaceAlpha',0.6)
% scatter(sparsemat(1:2:end,ii),...
%          deltamat(1:2:end,ii),30,hex2rgb(hexcolor{ii}),...
%          'filled','MarkerFaceAlpha',0.6)
hold on
end
hold off

zlim([0,350])
xlim([0.1,4])
ylim([0,0.5])
set(gca,'LineWidth',1,'FontSize',14,'Xscale','log');
box off
set(gcf,'color','w','PaperPositionMode', 'auto');
grid off
% saveas(gca, "panel_1E.svg")

end