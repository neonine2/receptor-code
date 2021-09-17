%% Panel 1A-1C
clear
% centers generated using DefineCellCenters(params.nCircles,params.rad);
load('cell_centers','centers');
% load ecm fiber network generated from combineECM.m
load('ecm','ecmMAT'); 

envparams = struct('rad',10,'nCircles',[3,8],'centers',centers,'Diff',45,...
              'dt',0.0178,'tTot',60,'trelease',2,'releaseC',7,...
              'csolUp',0.01,'mu',2.5,'hg',2,'fspeed',2,'gammas',0.001,...
              'gammab',0.00001,'kecm',9.3e-5,'ecmMAT',ecmMAT,'xecmpos',0,...
              'ecmsiteC',520,'uptakeD',1.15);
save('tissue_300by900_param','envparams')

%generate tissue env
sim_tissue(envparams,'tissue_300by900');

%generate soil env
% **run file named env_LGCP.R to generate the data env file**;

%% optimization
cellrad = [5,10,20];
MI_opt("tissue","tissue_300by900",cellrad); %tissue_300by900_szopt
MI_opt("grad","tissue_300by900",cellrad); %tissue_300by900_grad_szopt
MI_opt("soil","soil_var_2",cellrad,"optgrid",[20,60]); %soil_var_2_szopt

%plotting panel 1A
fnamelist = ["soil_var_2_szopt","tissue_300by900_szopt",...
                "tissue_300by900_grad_szopt"];
make_panel_1A(fnamelist)

%plotting panel 1B
fnamelist = ["soil_var_2_szopt","tissue_300by900_szopt",...
                           "tissue_300by900_grad_szopt"];
make_panel_1B(fnamelist)

%plotting panel 1C
fnamelist = ["soil_var_2_szopt","tissue_300by900_szopt",...
                          "tissue_300by900_grad_szopt"];
make_panel_1C(fnamelist)
                   
%% panel 1D
%generate tissue_kecm_1 to tissue_kecm_9
clear
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

%running optimization over all tissue env and fitted gradient
filelist = ["tissue_kecm_1","tissue_kecm_2","tissue_kecm_3",...
            "tissue_kecm_4","tissue_kecm_5","tissue_kecm_6",...
            "tissue_kecm_7","tissue_kecm_8","tissue_kecm_9"
            ]; 
nkecm = length(filelist);
cellrad = 10;

rel_eff_tissue = zeros(nkecm,1);
for ii = 1:nkecm
    %output file saved to current working directory
    rel_eff_tissue(ii) = MI_opt("tissue",filelist(ii),cellrad);
end
save('rel_eff_tissue_kcm','rel_eff_tissue','filelist')

rel_eff_grad = zeros(nkecm,1);
for ii = 1:nkecm
    rel_eff_grad(ii) = MI_opt("grad",filelist(ii),cellrad);
end
save('rel_eff_grad_kcm','rel_eff_grad','filelist')

%generate soil_var_xx
% **run file named soil_ge_var.R to generate the data env files**;

filelist = ["soil_var_0.001","soil_var_0.003",...
            "soil_var_0.01","soil_var_0.032",...
            "soil_var_0.1","soil_var_0.316",...
            "soil_var_1","soil_var_3.162",...
            "soil_var_10","soil_var_31.623","soil_var_80"]; 
            %var from logspace(-3,1.5,10) and 80;
nvar = length(filelist);

% running optimization over all soil env
cellrad = 10;
rel_eff_soil = zeros(nvar,1);
for ii = 1:nvar
    %output file saved to current working directory
    rel_eff_soil(ii) = ...
        MI_opt("soil",strcat(filelist(ii),".mat"),cellrad,"optgrid",[20,60]); 
end

save('soil_var_rel_eff','rel_eff_soil','filelist')

%plotting panel 1D
fnamelist = ["rel_eff_tissue_kcm","rel_eff_grad_kcm",...
                          "soil_var_rel_eff"];
make_panel_1D(fnamelist)
                            
%plotting panel 1E
clear
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
make_panel_1E(fnamelist)





