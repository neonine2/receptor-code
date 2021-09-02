cellrad = [5,10,20];
gridsz = [30,90];
MI_opt("tissue","tissue_300by900",cellrad,"optgrid",gridsz);
MI_opt("grad","tissue_300by900",cellrad,"optgrid",gridsz);
MI_opt("soil","soil_var_2",cellrad,"optgrid",gridsz);

load('colorpalette')
fnamelist = ["soil_var_2_szopt_large","tissue_300by900_szopt_large",...
                "tissue_300by900_grad_szopt_large"];
nenv = length(fnamelist);
load(fnamelist(2),'envmean','optMI','unifMI')
[nloc,m,ncell] = size(envmean);
sparsemat = zeros(nloc*ncell,nenv);
avgmat = zeros(nloc*ncell,nenv);
deltamat = zeros(nloc*ncell,nenv);
maxsparsity = vecnorm(ones(1,m),1,2)./vecnorm(ones(1,m),2,2);
% maxsparsity = -log(1/m);

for ii = 1:nenv
    load(fnamelist(ii),'envmean','optMI','unifMI')
    deltaI = (optMI-unifMI)./unifMI .* 100;
    %remove -ve values
    deltaI(deltaI<0) = 0/0; % make negative nan
    sparsity = 1-squeeze(vecnorm(envmean,1,2)./vecnorm(envmean,2,2))/maxsparsity;
%     sparsity = 1-sum(-envmean./sum(envmean,2).*log(envmean./sum(envmean,2)),2)/maxsparsity;
    avgconc = squeeze(mean(envmean,2));
    deltamat(:,ii) = deltaI(:);
    avgmat(:,ii) = avgconc(:);
    sparsemat(:,ii) = sparsity(:);
end

% making scatterplot
for jj = 1:nenv
ii = 4-jj;
jitterAmount = 0.3;
jitterValuesX = 2*(rand(size(avgmat(:,ii)))-0.5)*jitterAmount;   % +/-jitterAmount max
jitterValuesY = 2*(rand(size(deltamat(:,ii)))-0.5)*jitterAmount;   % +/-jitterAmount max
scatter3(avgmat(:,ii)+jitterValuesX,sparsemat(:,ii),...
         deltamat(:,ii),32,hex2rgb(hexcolor{ii}),...
         'filled','MarkerFaceAlpha',0.6)
        hold on
end
hold off

zlim([0,300])
xlim([0.05,2])
ylim([0,0.5])
set(gca,'LineWidth',1,'FontSize',14,'Xscale','log');
box off
set(gcf,'color','w','PaperPositionMode', 'auto');
grid off
saveas(gca, "panel_1E.svg")                