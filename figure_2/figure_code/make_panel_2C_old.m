function [] = make_panel_2C_old(fnamelist,mode)
%% panel C

X = categorical({'Soil','Tissue','Simple gradient'});
X = reordercats(X,{'Soil','Tissue','Simple gradient'});
load(fnamelist{1},'radlist')
nfile = length(fnamelist);
nrad = length(radlist);
deltaIMAT = zeros(nrad,nfile);
for ii = 1:nfile
    load(fnamelist{ii},'optMI','unifMI')
    deltaIMAT(:,ii) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
end

load('colorpalette')
b = bar(X,deltaIMAT','FaceColor','flat','EdgeColor','k','LineWidth',1.5);
for ii = 1:nfile
    for jj = 1:nrad
        b(jj).CData(ii,:) = hex2rgb(hexcolor{ii});
    end
end
set(gca,'Linewidth',1.5,'xtick',[],'fontsize',16)
box off
pbaspect([3.2,1,1])
ylabel('$\eta$','Interpreter','latex','fontsize',21)
saveas(gca,"panel_1C.png")

end