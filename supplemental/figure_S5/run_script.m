%panel B
fnamelist = ["soil_var_2_szopt","tissue_300by900_szopt",...
                "tissue_300by900_grad_szopt"];
make_panel_B(fnamelist,"ratio")
make_panel_B(fnamelist,"diff")
make_panel_B(fnamelist,"unif")

%panel C
tissue_filelist = ["tissue_kecm_1","tissue_kecm_2","tissue_kecm_3",...
                    "tissue_kecm_4","tissue_kecm_5","tissue_kecm_6",...
                    "tissue_kecm_7","tissue_kecm_8","tissue_kecm_9"]; 
        
soil_filelist = ["soil_var_0.001","soil_var_0.003",...
                    "soil_var_0.01","soil_var_0.032",...
                    "soil_var_0.1","soil_var_0.316",...
                    "soil_var_1","soil_var_3.162",...
                    "soil_var_10","soil_var_31.623","soil_var_80"]; 

make_panel_2D(tissue_filelist,soil_filelist,"ratio")
make_panel_2D(tissue_filelist,soil_filelist,"diff")
make_panel_2D(tissue_filelist,soil_filelist,"unif")

%panel D
fnamelist = ["soil_var_2_szopt_large",...
                "tissue_300by900_szopt_large",...
                "tissue_300by900_grad_szopt_large"];
% plotting panel 2E
make_panel_2E(fnamelist,"ratio")
make_panel_2E(fnamelist,"diff")
make_panel_2E(fnamelist,"unif")