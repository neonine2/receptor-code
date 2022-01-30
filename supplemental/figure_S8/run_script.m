clear all;
close all;

fname = ["tissue_300by900_szopt.mat","soil_var_2_szopt.mat"];
testing_param = ["koff","h"];
mode = ["static","dynamic"];
cellmodel = 1; %5um
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

%% plotting

% IMPORTANT: make sure to rename file as this will generate panel5 name
% file
make_panel_5C(["tissue_300by900_szopt_scheme_h_static_5um",...
                "soil_var_2_szopt_scheme_h_static_5um",...
                  "tissue_300by900_szopt_scheme_koff_static_5um",...
                     "soil_var_2_szopt_scheme_koff_static_5um"])
                 
make_panel_5F(["tissue_300by900_szopt_scheme_h_dynamic_5um",...
                     "soil_var_2_szopt_scheme_h_dynamic_5um",...
                     "tissue_300by900_szopt_scheme_koff_dynamic_5um",...
                     "soil_var_2_szopt_scheme_koff_dynamic_5um"])