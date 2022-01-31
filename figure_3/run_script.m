clear all;
close all;
mode = "diff";
% perturb optimal receptor placement and compute efficiency in tissue
cellmodel = 2; % 10um radius
fname = which("tissue_300by900_szopt.mat");
assess_perturb(fname,cellmodel) % tissue_300by900_szopt_perturb_10um
% plotting panel B
make_panel_3B(fname,which("tissue_300by900_szopt_perturb_10um.mat"),'mode',mode)

% perturb optimal receptor placement and compute efficiency in soil
fname = which("soil_var_2_szopt.mat");
assess_perturb(fname,cellmodel) % soil_var_2_szopt_perturb_10um
% plotting panel B
make_panel_3B(fname,which("soil_var_2_szopt_perturb_10um.mat"),'mode',mode)

