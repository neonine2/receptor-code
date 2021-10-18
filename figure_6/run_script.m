%% making navigation environment for localization task
clear
load('tissue_300by900_param','envparams');
envparams.nCircles = [1,2]; %consider a smaller arena, note these cells are
                            %simply placed to obtained an appropriately
                            %sized arena, they do not interact with the
                            %moving cell
envparams.centers = [20,-400;20,400];
envparams.tTot = 100;
envparams.xecmpos = 30; %changing this parameter sets the part of the ecm
                        %network to use in generating the tissue env

% get first env
fname = "tissue_taxis_1";
sim_tissue(envparams,fname)

% get second env
fname = "tissue_taxis_2";
envparams.xecmpos = 50;
sim_tissue(envparams,fname) % simulate env

%% setting paramaters
clear
load('tissue_300by900_szopt','receptor_params')
load('scheme_parameter')

% setting feedback scheme parameter
thour = 6;
scheme_param.mean_cell_radius = 10;
conversion_factor = conc2count(scheme_param.mean_cell_radius,scheme_param.N);
scheme_param.T = 60*60*thour; % seconds
scheme_param.kd = receptor_params.kd*conversion_factor;
scheme_param.rtot = receptor_params.rtot;
scheme_param.receptornoise = receptor_params.receptornoise;

% default parameter
% let cell navigate over the default tissue env as well as the two
% generated above, main difference between these env are ecm network
% pattern
fnamelist = ["tissue_300by900","tissue_taxis_1","tissue_taxis_2"];
for ii = 1:length(fnamelist)
    fname = fnamelist(ii);
    racing_cells(fname, scheme_param, "envmodel","tissue",...
        "task","localization");
    racing_cells(fname, scheme_param, "envmodel","grad",...
        "task","localization");
end

%% making tissue environment for retention task
clear all;
close all;

load('tissue_300by900_param','envparams');
envparams.nCircles = [1,2];
envparams.centers = [20,-400;20,400];
envparams.Diff = 5;
envparams.releaseC = 20;
envparams.fspeed = 0.3;
envparams.gammas = 0.3;
envparams.gammab = 0.003;
envparams.dt = 0.16;
envparams.tTot = 40;

% adjust xecmpos to get different ecm pattern;
envparams.xecmpos = 20; 
fname = "tissue_retention_1";
sim_tissue(envparams,fname)

envparams.xecmpos = 40;
fname = "tissue_retention_2";
sim_tissue(envparams,fname)

envparams.xecmpos = 60;
fname = "tissue_retention_3";
sim_tissue(envparams,fname)

%% setting paramaters
clear
load('tissue_300by900_szopt','receptor_params')
load('scheme_parameter')

% setting feedback scheme parameter
thour = 2;
scheme_param.mean_cell_radius = ellipsePerimeter(2,5)/(2*pi);
conversion_factor = conc2count(scheme_param.mean_cell_radius,scheme_param.N);
scheme_param.T = 60*60*thour; % seconds
scheme_param.kd = receptor_params.kd*conversion_factor;
scheme_param.rtot = receptor_params.rtot;
scheme_param.receptornoise = receptor_params.receptornoise;

% default parameter
fnamelist = ["tissue_retention_1","tissue_retention_2","tissue_retention_3"];
for ii = 1:length(fnamelist)
    fname = fnamelist(ii);
    racing_cells(fname,scheme_param,"envmodel","tissue",...
        "task","retention");
    racing_cells(fname,scheme_param,"envmodel","grad",...
        "task","retention");
end

%% plotting
%localization task
make_panel_6B("tissue_300by900_szopt.mat");
make_panel_6C("tissue_300by900_localization_feedback");

fnamelist = ["tissue_300by900_localization_feedback",...
    "tissue_taxis_1_localization_feedback",...
    "tissue_taxis_2_localization_feedback"];
fnamelist_grad = ["tissue_300by900_grad_localization_feedback",...
    "tissue_taxis_1_grad_localization_feedback",...
    "tissue_taxis_2_grad_localization_feedback"];
make_panel_6D_6E(fnamelist,fnamelist_grad);

%retention task
make_panel_6G("tissue_300by900_szopt.mat");
make_panel_6H("tissue_retention_1_retention_feedback");

fnamelist = ["tissue_retention_1_retention_feedback",...
                "tissue_retention_2_retention_feedback",...
                "tissue_retention_3_retention_feedback"];
fnamelist_grad = ["tissue_retention_1_grad_retention_feedback",...
                    "tissue_retention_2_grad_retention_feedback",...
                    "tissue_retention_3_grad_retention_feedback"];
make_panel_6I_6J(fnamelist,fnamelist_grad);

