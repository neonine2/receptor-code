cellmodel = 2; % 10um radius
fname = which("tissue_300by900_szopt.mat");
assess_perturb(fname,cellmodel)
% plotting panel B
make_panel_3B(fname,which("tissue_300by900_szopt_perturb_10um.mat"))

fname = which("soil_var_2_szopt.mat");
assess_perturb(fname,cellmodel)
% plotting panel B
make_panel_3B(fname,which("soil_var_2_szopt_perturb_10um.mat"))

