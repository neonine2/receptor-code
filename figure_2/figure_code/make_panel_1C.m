function [] = make_panel_1C(fnamelist)
%% panel C
X = categorical({'Soil','Tissue','Simple gradient'});
X = reordercats(X,{'Soil','Tissue','Simple gradient'});
load(fnamelist{1},'radlist')
nfile = length(fnamelist);
nrad = length(radlist);
deltaIMAT = zeros(nrad,nfile);
for ii = 1:nfile
    load(fnamelist{ii},'optMI','unifMI')
    deltaIMAT(:,ii) = (sum(optMI) - sum(unifMI))./sum(unifMI) .* 100;
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
ylabel('$\Delta I/I$','Interpreter','latex','fontsize',21)

end