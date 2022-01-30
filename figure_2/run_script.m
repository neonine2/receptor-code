clear all;
close all;

%% generating environments (ligand field)
envparams = struct('rad',10,'nCircles',[3,8],'Diff',45,...
              'dt',0.0178,'tTot',60,'trelease',2,'releaseC',7,...
              'csolUp',0.01,'mu',2.5,'hg',2,'fspeed',2,'gammas',0.001,...
              'gammab',0.00001,'kecm',9.3e-5,'xecmpos',0,...
              'ecmsiteC',520,'uptakeD',1.15);
          
%generating cells in environment
% cell positions are stochastically generated, so expect some variability
% centers = generateCellPos(envparams.nCircles,envparams.rad); 
% save('cell_centers','centers');
load('cell_centers','centers');
envparams.centers = centers;

% generating ecm fiber network
% ecm network is stochastically generated, so expect some variability
% ecm = generateECM(40000);
% [ecmMAT,~,~] = histcounts2(ecm(:,1),ecm(:,2),100:2300,100:1100);
% imagesc(ecmMAT')
% colorbar()
% save('ecm.mat','ecmMAT','-v7.3')
load('ecm','ecmMAT'); 
envparams.ecmMAT = ecmMAT;

%saving all environmental parameters
save('tissue_300by900_param','envparams')

%generate default tissue env
sim_tissue(envparams,'tissue_300by900');

%generate default soil env
% IMPORTANT: run file named env_LGCP.R to generate the data env file, this
% will again be stochastic, expect variability

%% panel 2A
cellrad = [5,10,20];
MI_opt("tissue","tissue_300by900",cellrad); %tissue_300by900_szopt
MI_opt("grad","tissue_300by900",cellrad); %tissue_300by900_grad_szopt
MI_opt("soil","soil_var_2",cellrad,"optgrid",[20,60]); %soil_var_2_szopt

%output files
fnamelist = ["soil_var_2_szopt","tissue_300by900_szopt",...
                "tissue_300by900_grad_szopt"];

%plotting panel 2A
make_panel_2A(fnamelist)

%plotting panel 2B
make_panel_2B(fnamelist)


%% panel 2C

%% running optimization
alpha_run("diff");
%plotting panel 2C
make_panel_2C("opt_result_diff.mat"); 
                   
%% panel 2D
%generate tissue with different kecm parameter
clear;
load('tissue_300by900_param','envparams');
nkecm = 9;
kecmlist = logspace(-4,-2,nkecm)*12/520;
envparams.trelease = 100;
envparams.releaseC = 0.14;
parfor ii = 1:nkecm
    par = envparams;
    par.kecm = kecmlist(ii);
    sim_tissue(par,strcat('tissue_kecm_',num2str(ii)))
end
% output tissue env
tissue_filelist = ["tissue_kecm_1","tissue_kecm_2","tissue_kecm_3",...
            "tissue_kecm_4","tissue_kecm_5","tissue_kecm_6",...
            "tissue_kecm_7","tissue_kecm_8","tissue_kecm_9"
            ]; 
cellrad = 10;

% optimize in tissue
rel_eff_tissue = zeros(nkecm,1);
for ii = 1:nkecm
    %output file saved to current working directory
    rel_eff_tissue(ii) = MI_opt("tissue",tissue_filelist(ii),cellrad);
end
save('rel_eff_tissue_kcm','rel_eff_tissue','filelist')

% optimize in fitted gradient
rel_eff_grad = zeros(nkecm,1);
for ii = 1:nkecm
    rel_eff_grad(ii) = MI_opt("grad",tissue_filelist(ii),cellrad);
end
save('rel_eff_grad_kcm','rel_eff_grad','filelist')

%generate soil with different ligand spatial spread
% **run file named soil_ge_var.R to generate the data env files**;
soil_filelist = ["soil_var_0.001","soil_var_0.003",...
            "soil_var_0.01","soil_var_0.032",...
            "soil_var_0.1","soil_var_0.316",...
            "soil_var_1","soil_var_3.162",...
            "soil_var_10","soil_var_31.623","soil_var_80"]; 
            %var from logspace(-3,1.5,10) and 80;
nvar = length(soil_filelist);
% running optimization over all soil env
cellrad = 10;
rel_eff_soil = zeros(nvar,1);
for ii = 1:nvar
    %output file saved to current working directory
    rel_eff_soil(ii) = ...
        MI_opt("soil",strcat(soil_filelist(ii),".mat"),cellrad,"optgrid",[20,60]); 
end
save('soil_var_rel_eff','rel_eff_soil','filelist')

%plotting panel 2D
make_panel_2D(tissue_filelist,soil_filelist,"diff")

%% panel 2E
clear
% optimize over finer grid to obtain more sample points
cellrad = [5,10,20];
gridsz = [30,90];
MI_opt("tissue","tissue_300by900",cellrad,"optgrid",gridsz,...
    "savefname","tissue_300by900_szopt_large");
MI_opt("grad","tissue_300by900",cellrad,"optgrid",gridsz,...
    "savefname","tissue_300by900_grad_szopt_large");
MI_opt("soil","soil_var_2",cellrad,"optgrid",gridsz,...
    "savefname","soil_var_2_szopt_large");
fnamelist = ["soil_var_2_szopt_large",...
                "tissue_300by900_szopt_large",...
                "tissue_300by900_grad_szopt_large"];
% plotting panel 2E
make_panel_2E(fnamelist,"diff")

