function [] = alpha_run(mode)

cellrad = [5,10];
nalpha = 10;
alphavec = logspace(-3,0,nalpha);
miMAT = zeros(length(alphavec),length(cellrad));
migradMAT = zeros(length(alphavec),length(cellrad));
misoilMAT = zeros(length(alphavec),length(cellrad));
for ii = 1:10
    receptor_params = struct('rtot',1000,'kd',40,...
                                'receptornoise',alphavec(ii));
    fname = which(strcat("alpha/tissue_300by900_szopt",num2str(ii),".mat"));
    disp(fname);
    if isfile(fname)
        load(fname,"optMI","unifMI")
        miMAT(ii,:) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    else
        miMAT(ii,:) = MI_opt("tissue","tissue_300by900",cellrad,...
            "receptor_params",receptor_params,"index",ii,"mode",mode); 
    end
    
    fname = which(strcat("alpha/tissue_300by900_grad_szopt",num2str(ii),".mat"));
    if isfile(fname)
        load(fname,"optMI","unifMI")
        migradMAT(ii,:) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    else
        migradMAT(ii,:) = MI_opt("grad","tissue_300by900",cellrad,...
             "receptor_params",receptor_params,"index",ii,"mode",mode);
    end
    
    fname = which(strcat("alpha/soil_var_2_szopt",num2str(ii),".mat"));
    if isfile(fname)
        load(fname,"optMI","unifMI")
        misoilMAT(ii,:) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    else
        misoilMAT(ii,:) = MI_opt("soil","soil_var_2",cellrad,...
             "receptor_params",receptor_params,"optgrid",[20,60],"index",ii,"mode",mode);
    end
end
save(strcat("opt_result_",mode,".mat"))
end

