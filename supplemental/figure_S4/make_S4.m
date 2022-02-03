function [] = make_S4(mode)
load(strcat("opt_result_",mode,".mat"))
load("colorpalette")
tiledlayout(1,2)
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
    set(gca,"fontsize",16,"Xscale","log")
    xlim([10^-2.5, 10^-0.4])
%     ylim([0.01,2])
    pbaspect([1,1,1])
    xlabel("$\alpha$",'interpreter','latex','fontsize',20)
    ylabel("$\eta$",'interpreter','latex','fontsize',20)
    title(strcat("Cell radius = ",num2str(cellrad(ii)),"\mu m"),...
        "FontWeight",'normal')
end
saveas(gca,strcat("fig_S4_",mode,".svg"))
end

