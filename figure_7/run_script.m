clear all;
close all;

%% running optimization for different receptor parameters
% this run was divided into two parts in order to speed up the process

% part 1
rtotlist = logspace(2,5,16);
kdlist = logspace(1,5,21);
paramlist = combvec(rtotlist,kdlist)';
nparam = length(paramlist);
rel_eff_mat = zeros(nparam,1);
cellrad = 10;
default_param = struct('rtot',1000,'kd',40,'receptornoise',0.1);
for jj = 1:nparam
    % setting parameters
    newparam = default_param;
    newparam.rtot = paramlist(jj,1);
    newparam.kd = paramlist(jj,2);
    rel_eff_mat(jj) = MI_opt("tissue","tissue_300by900",cellrad,...
        "receptor_params",newparam,"saving_data",true,"plot_env",false,...
        "index",jj);
    if mod(jj, 10) == 0
        disp(jj);
    end
end 
save("kd_rtot_rel_eff");

% part 2 runs over parameters where the optimization is very slow, so we
% use a coarser grid
clear
rtotlist = logspace(2,5,7);
kdlist = logspace(-2,0.5,6);
paramlist = combvec(rtotlist,kdlist)';
nparam = length(paramlist);
rel_eff_mat = zeros(nparam,1);
cellrad = 10;
default_param = struct('rtot',1000,'kd',40,'receptornoise',0.1);
for jj = 1:nparam
    % setting parameters
    newparam = default_param;
    newparam.rtot = paramlist(jj,1);
    newparam.kd = paramlist(jj,2);
    rel_eff_mat(jj) = MI_opt("tissue","tissue_300by900",cellrad,...
        "receptor_params",newparam,"saving_data",true,"plot_env",false,...
        "index",jj);
    if mod(jj, 10) == 0
        disp(jj);
    end
end 
save("kd_rtot_rel_eff_2");

%% plotting
fnamelist = ["kd_rtot_rel_eff","kd_rtot_rel_eff_2"];
make_panel_7B(fnamelist)

