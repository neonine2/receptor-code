clear all;
close all;

% testing efficiency of feedback scheme in both static and dynamic
% environment and across different parameters
fname = ["tissue_300by900_szopt.mat","soil_var_2_szopt.mat"]; %tissue and soil
testing_param = ["koff","h"]; % endocytosis and incorporation rate
mode = ["static","dynamic"]; % static and dynamic cell, in the dynamic case
                             % the cell is translated along a straight path
                             % at a speed of 2um/min
cellmodel = 2; %10um radius cell
% load parameters for feedback scheme simulation
load('scheme_parameter','scheme_param') 

for ii = 1:length(fname)
    for jj = 1:length(testing_param)
        for kk = 1:length(mode)
            feedback_scheme_efficiency(which(fname(ii)),scheme_param,...
                                    testing_param(jj),mode(kk),cellmodel)
        end
    end
end

%% plotting 5B, 5C, 5F
make_panel_5B(["tissue_300by900_szopt_scheme_koff_static_10um",...
                "soil_var_2_szopt_scheme_koff_static_10um"])

make_panel_5C(["tissue_300by900_szopt_scheme_h_static_10um",...
                "soil_var_2_szopt_scheme_h_static_10um",...
                  "tissue_300by900_szopt_scheme_koff_static_10um",...
                     "soil_var_2_szopt_scheme_koff_static_10um"])
                 
make_panel_5F(["tissue_300by900_szopt_scheme_h_dynamic_10um",...
                     "soil_var_2_szopt_scheme_h_dynamic_10um",...
                     "tissue_300by900_szopt_scheme_koff_dynamic_10um",...
                     "soil_var_2_szopt_scheme_koff_dynamic_10um"])
                 
%% generate data for kymograph
clear

% translate cell in a straight line down the y direction for nminute at a
% speed of 2um/minute
load('scheme_parameter')
load('tissue_300by900_szopt','centerlist','fconc','fenv',...
                'conversion_factor','receptor_params','radlist');
cellmodel = 2;
scheme_param.mean_cell_radius = radlist(cellmodel);
fcount = @(x,y) fconc(x,y)*conversion_factor(cellmodel);
scheme_param.fcount = fcount;
nminute = 360; % 6 hours
scheme_param.T = 60*nminute; 
scheme_param.kd = receptor_params.kd*conversion_factor(cellmodel);
scheme_param.rtot = receptor_params.rtot;
scheme_param.receptornoise = receptor_params.receptornoise;
scheme_param.h = 0.004;
scheme_param.koff = 0.18;

% set cell shape and starting position
m = scheme_param.N;
phi = linspace(pi,-pi+(2*pi)/m,m)';
phi = circshift(flipud(phi),m/2+1);
phi(1) = 0; % correct for minor numerical inaccuracy
s = scheme_param.mean_cell_radius;
scheme_param.cellboundary = s*[cos(phi),sin(phi)];
scheme_param.dx = ellipsePerimeter(s,s)/m;
scheme_param.cellp = centerlist(2,:,cellmodel);

% move cell and obtain receptor profiles generated due to feedback scheme
results = cell_move(scheme_param,"mode","dynamic",...
                 "decoder_method","yaxis1d",'receptor',"feedback");
save('kymograph_data')
disp(mean(results.stat,1)') %show some output data from feedback scheme

%% plotting 5D, 5E and 5G
make_panel_5D_5E_5G("kymograph_data")
    