function [] = make_panel_2C(fname)

load(fname)
load("colorpalette")
tiledlayout(1,2)
alphavec = logspace(-3,0,nalpha);
for ii = 1:2
    nexttile
    plot(alphavec, misoilMAT(:,ii),"color",hexcolor{1},...
        "linewidth",1.5);
    hold on
    plot(alphavec, miMAT(:,ii),"color",hexcolor{2},...
        "linewidth",1.5);
    plot(alphavec, migradMAT(:,ii),"color",hexcolor{3},...
        "linewidth",1.5);
    hold off
    set(gca,"fontsize",16,"Xscale","log","Yscale","log")
    xlim([10^-2.5, 10^-0.6]) %10^-0.6
    ylim([0.01,2])
    pbaspect([2,1,1])
    xlabel("$\alpha$",'interpreter','latex','fontsize',20)
    ylabel("$\eta$",'interpreter','latex','fontsize',20)
    title(strcat("Cell radius = ",num2str(cellrad(ii)),"\mu m"),...
        "FontWeight",'normal')
end
saveas(gca,"fig_2C.svg")

end

