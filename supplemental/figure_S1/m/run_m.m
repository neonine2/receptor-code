%% running optimization
clear all;
close all;

mode = "diff";
cellrad = [5,10];
nbin = 6;
binvec = linspace(50,300,nbin);
miMAT = zeros(length(binvec),length(cellrad));
migradMAT = zeros(length(binvec),length(cellrad));
misoilMAT = zeros(length(binvec),length(cellrad));
for ii = 1:nbin
    receptor_params = struct('rtot',1000,'kd',40,'receptornoise',0.01);
    fname = strcat("tissue_300by900_szopt",num2str(ii),".mat");
    if isfile(fname)
        load(fname,"optMI","unifMI")
        miMAT(ii,:) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    else
        miMAT(ii,:) = MI_opt("tissue","tissue_300by900",cellrad,...
            "index",ii,"mode",mode,"nbin",binvec(ii),"receptor_params",receptor_params); 
    end
    
    fname = strcat("tissue_300by900_grad_szopt",num2str(ii),".mat");
    if isfile(fname)
        load(fname,"optMI","unifMI")
        migradMAT(ii,:) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    else
        migradMAT(ii,:) = MI_opt("grad","tissue_300by900",cellrad,...
            "index",ii,"mode",mode,"nbin",binvec(ii),"receptor_params",receptor_params);
    end
    
    fname = strcat("soil_var_2_szopt",num2str(ii),".mat");
    if isfile(fname)
        load(fname,"optMI","unifMI")
        misoilMAT(ii,:) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    else
        misoilMAT(ii,:) = MI_opt("soil","soil_var_2",cellrad,...
                "optgrid",[20,60],"index",ii,"mode",mode,"nbin",binvec(ii),...
                                        "receptor_params",receptor_params);
    end
end
save(strcat("opt_result_",mode,".mat"))

%% plotting Fig2C
clear
load("opt_result_diff.mat")
load("colorpalette")
tiledlayout(1,2)
for ii = 1:2
    nexttile
    plot(binvec, misoilMAT(:,ii),"color",hexcolor{1},...
        "linewidth",1.5);
    hold on
    plot(binvec, miMAT(:,ii),"color",hexcolor{2},...
        "linewidth",1.5);
    plot(binvec, migradMAT(:,ii),"color",hexcolor{3},...
        "linewidth",1.5);
    hold off
    set(gca,"fontsize",12,"Xscale","log")
    xlim([100, 300]) %10^-0.6
%     ylim([0.01,2])
    pbaspect([1,1,1])
    xlabel("m",'fontsize',14)
    ylabel("$\eta$",'interpreter','latex','fontsize',16)
    title(strcat("Cell radius = ",num2str(cellrad(ii)),"\mu m"),...
        "FontWeight",'normal')
    legend('Soil','Tissue','Gradient','location','northwest','box','off')
end
saveas(gca,"membin_efficacy.svg")
