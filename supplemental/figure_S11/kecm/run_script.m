clear all
close all

%% compute efficacy for real receptors across different alpha

%list of real receptor parameters
receptornames = {'IL-2R','TNFR-1','TGFbR-2','CR3','C5aR','CCR2','CXCR4',...
                                                  'CCR5','GABAAR','Robo1'};
rtot = [1000,1000,10000,40000,50000,1100,1572,593,200,22300];
kd = [0.013,0.019,0.05,12.5,2,1.53,5,4,12,235];
nreceptor = 10;

%list of kcem (tissue env) to scan over
envlist = ["tissue_kecm_1","tissue_kecm_2","tissue_kecm_3",...
            "tissue_kecm_4","tissue_kecm_5"];
nenv = length(envlist);

%computing optimal efficacy
rel_eff_mat = zeros(nreceptor,nenv);
cellrad = 10;
mode = "ratio";
default_param = struct('rtot',1000,'kd',40,'receptornoise',0.1);
for jj = 1:nreceptor
    for ii = 1:nenv
        fname = strcat(envlist(ii),"_opt",num2str((ii-1)*nreceptor+jj),".mat");
        if isfile(fname)
            load(fname,"optMI","unifMI")
            rel_eff_mat(jj,ii) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
            continue
        end
        % setting parameters
        newparam = default_param;
        newparam.rtot = rtot(jj);
        newparam.kd = kd(jj);
        rel_eff_mat(jj,ii) = MI_opt("tissue",envlist(ii),cellrad,...
            "receptor_params",newparam,"saving_data",true,"plot_env",false,...
            "index",(ii-1)*nreceptor+jj,"cfac",newparam.kd/20,"optgrid",[5,15]);
    end 
    disp(jj);
end

%plotting results
tiledlayout(nenv,1)
for ii = 1:nenv
    nexttile
    b = bar(rel_eff_mat(:,ii));
    b.FaceColor = 'flat';
    b.CData(6:end,:) = repmat([1 0 0],5,1);
    ylabel("\fontsize{6} \eta")
    set(gca,"fontsize",6)
    pbaspect([2,1,1])
end
set(gca,'xticklabel',receptornames)
xtickangle(45)
saveas(gca,strcat("kecm_",mode),'png')