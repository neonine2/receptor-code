clear

fname = ["tissue_300by900_szopt.mat","soil_var_2_szopt.mat"];
testing_param = ["koff","h"];
mode = ["static","dynamic"];
cellmodel = 1; %5um

for ii = 1:2
    for jj = 1:2
        for kk = 1:2
            feedback_scheme_efficiency(which(fname(ii)),testing_param(jj),mode(kk),cellmodel)
        end
    end
end

%% plotting

make_panel_4C(["tissue_300by900_szopt_scheme_h_static_5um",...
                "soil_var_2_szopt_scheme_h_static_5um",...
                  "tissue_300by900_szopt_scheme_koff_static_5um",...
                     "soil_var_2_szopt_scheme_koff_static_5um"])
                 
make_panel_4G(["tissue_300by900_szopt_scheme_h_dynamic_5um",...
                     "soil_var_2_szopt_scheme_h_dynamic_5um",...
                     "tissue_300by900_szopt_scheme_koff_dynamic_5um",...
                     "soil_var_2_szopt_scheme_koff_dynamic_5um"])