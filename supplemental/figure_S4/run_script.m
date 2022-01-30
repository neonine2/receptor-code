%% generate data files
alpha_run("diff");
alpha_run("unif");
alpha_run("expdiff");
alpha_run("ratio");

%% plotting FigS4
make_S4("diff");
make_S4("unif");
make_S4("expdiff");
make_S4("ratio");

tissue_filelist = ["tissue_kecm_1","tissue_kecm_2","tissue_kecm_3",...
            "tissue_kecm_4","tissue_kecm_5","tissue_kecm_6",...
            "tissue_kecm_7","tissue_kecm_8","tissue_kecm_9"
            ]; 
soil_filelist = ["soil_var_0.001","soil_var_0.003",...
            "soil_var_0.01","soil_var_0.032",...
            "soil_var_0.1","soil_var_0.316",...
            "soil_var_1","soil_var_3.162",...
            "soil_var_10","soil_var_31.623","soil_var_80"]; 
make_panel_2D(tissue_filelist,soil_filelist,"ratio")

fnamelist = ["soil_var_2_szopt_large",...
                "tissue_300by900_szopt_large",...
                "tissue_300by900_grad_szopt_large"];
% plotting panel 2E
make_panel_2E(fnamelist,"expdiff")