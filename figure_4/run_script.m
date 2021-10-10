cellmodel = 2; % 10um radius
fname = which("tissue_300by900_szopt.mat");
assess_perturb(fname,cellmodel)
% plotting panel B
make_panel_4B(fname)