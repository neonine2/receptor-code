%% running optimization
clear
cellrad = [5,10,20];
nalpha = 10;
alphavec = logspace(-3,-0.5,nalpha);

miMAT = zeros(nalpha,length(cellrad));
migradMAT = zeros(nalpha,length(cellrad));
misoilMAT = zeros(nalpha,length(cellrad));
load("opt_result.mat")
for ii = 1:nalpha
    receptor_params = struct('rtot',1000,'kd',40,...
                                'receptornoise',alphavec(ii));
    miMAT(ii,:) = MI_opt("tissue","tissue_300by900",cellrad,...
        "receptor_params",receptor_params,'index',ii); 
    migradMAT(ii,:) = MI_opt("grad","tissue_300by900",cellrad,...
        "receptor_params",receptor_params,'index',ii);
    misoilMAT(ii,:) = MI_opt("soil","soil_var_2",cellrad,...
        "receptor_params",receptor_params,'index',ii);
end
save("opt_result.mat")

%% plotting
clear
load("opt_result.mat")
load("colorpalette")
tiledlayout(1,2)
for ii = 1:2
    nexttile
    plot(alphavec, misoilMAT(:,ii),"color",hexcolor{1},...
        "linewidth",1);
    hold on
    plot(alphavec, miMAT(:,ii),"color",hexcolor{2},...
        "linewidth",1);
    plot(alphavec, migradMAT(:,ii),"color",hexcolor{3},...
        "linewidth",1);
    hold off
    set(gca,"fontsize",16,"Xscale","log","Yscale","log")
    xlim([10^-3, 10^-0.6])
    pbaspect([1,1,1])
    xlabel("$\alpha$",'interpreter','latex','fontsize',20)
    ylabel("$\Delta I/I$",'interpreter','latex','fontsize',20)
    title(strcat("Cell radius = ",num2str(cellrad(ii)),"\mu m"),...
        "FontWeight",'normal')
    legend('Soil','Tissue','Gradient','linewidth',1,...
        'Location','southeast','fontsize',14)
    legend('boxoff')
end
saveas(gca,"alpha_deltaI.svg")





