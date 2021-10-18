cellmodel = 1; % 5um radius
fname = which("tissue_300by900_szopt.mat");
assess_perturb(fname,cellmodel)
make_panel_B(fname,which("tissue_300by900_szopt_perturb_5um.mat"))

fname = which("soil_var_2_szopt.mat");
assess_perturb(fname,cellmodel)
make_panel_B(fname,which("soil_var_2_szopt_perturb_5um.mat"))