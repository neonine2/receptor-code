clear all
close all

%% compute efficacy for real receptors across different alpha

%list of real receptor parameters
receptornames = {'IL-2R','TNFR-1','TGF\beta R-2','CR3','C5aR','CCR2','CXCR4',...
                                                  'CCR5','GABA_AR','Robo1'};
rtot = [1000,1000,10000,40000,50000,1100,1572,593,200,22300];
kd = [0.013,0.019,0.05,12.5,2,1.53,5,4,12,235];
nreceptor = 10;

%list of alpha to scan over
nparam = 5;
alphalist = logspace(-3,-1,nparam);

%computing optimal efficacy
rel_eff_mat = zeros(nreceptor,nparam);
cellrad = 10;
mode = "diff";
default_param = struct('rtot',1000,'kd',40,'receptornoise',0.1);
for jj = 1:nreceptor
    for ii = 1:nparam
        fname = strcat("tissue_300by900_opt",num2str((ii-1)*nreceptor+jj),".mat");
        if isfile(fname)
            load(fname,"optMI","unifMI")
            rel_eff_mat(jj,ii) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
            continue
        end
        % setting parameters
        newparam = default_param;
        newparam.rtot = rtot(jj);
        newparam.kd = kd(jj);
        newparam.receptornoise = alphalist(ii);
        rel_eff_mat(jj,ii) = MI_opt("tissue","tissue_300by900",cellrad,...
            "receptor_params",newparam,"saving_data",true,"plot_env",false,...
            "index",(ii-1)*nreceptor+jj,"cfac",newparam.kd/20,"optgrid",[5,15]);
    end 
    disp(jj);
end

%plotting results
tiledlayout(nparam,1)
for ii = 1:nparam
    nexttile
    b = bar(rel_eff_mat(:,ii));
    b.FaceColor = 'flat';
    b.CData(6:end,:) = repmat([1 0 0],5,1);
    ylabel("\fontsize{11} \eta")
    set(gca,"fontsize",8)
    pbaspect([2,1,1])
end
set(gca,'xticklabel',receptornames,'fontsize',8,'linewidth',0.5)
xtickangle(45)
% saveas(gca,strcat("alpha_",mode),'svg')