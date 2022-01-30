clear all;
close all;

%% running optimization for different receptor parameters
% this run was divided into two parts in order to speed up the process

% part 1
mode = "ratio";
rtotlist = logspace(2,5,16);
kdlist = logspace(1,5,21);
kdlist = kdlist(1:9);
paramlist = combvec(rtotlist,kdlist)';
nparam = length(paramlist);
rel_eff_mat = zeros(nparam,1);
cellrad = 10;
default_param = struct('rtot',1000,'kd',40,'receptornoise',0.1);
for jj = 1:nparam
    if isfile(strcat("tissue_300by900_opt",num2str(jj),".mat"))
        continue
    end
    % setting parameters
    newparam = default_param;
    newparam.rtot = paramlist(jj,1);
    newparam.kd = paramlist(jj,2);
    rel_eff_mat(jj) = MI_opt("tissue","tissue_300by900",cellrad,...
        "receptor_params",newparam,"saving_data",true,"plot_env",false,...
        "index",jj,"cfac",newparam.kd/20,"optgrid",[5,15],"mode",mode);
    if mod(jj, 10) == 0
        disp(jj);
    end
end 
save("kd_rtot_rel_eff");

% part 2 runs over parameters where the optimization is very slow, so we
% use a coarser grid
clear
rtotlist = logspace(2,3.5,8);
kdlist = logspace(-2,-1,6);
paramlist = combvec(rtotlist,kdlist)';
nparam = length(paramlist);
rel_eff_mat = zeros(nparam,1);
cellrad = 10;
default_param = struct('rtot',1000,'kd',40,'receptornoise',0.1);
parfor jj = 1:nparam
    if isfile(strcat("tissue_300by900_opt",num2str(jj),".mat"))
        continue
    end
    % setting parameters
    newparam = default_param;
    newparam.rtot = paramlist(jj,1);
    newparam.kd = paramlist(jj,2);
    rel_eff_mat(jj) = MI_opt("tissue","tissue_300by900",cellrad,...
        "receptor_params",newparam,"saving_data",true,"plot_env",false,...
        "index",jj,"cfac",newparam.kd/20,"optgrid",[5,15],"mode",mode);
    if mod(jj, 10) == 0
        disp(jj);
    end
end 
save("kd_rtot_rel_eff_2");

%% plotting
clear all;
mode = "diff";
rel_eff_mat = NaN(144,1);
for ii = 1:length(rel_eff_mat)
    fname = strcat("part1/tissue_300by900_opt",num2str(ii),".mat");
    if exist(fname,'file')~=0
        load(fname,"optMI","unifMI")
        rel_eff_mat(ii) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    end
end
rtotlist = logspace(2,5,16);
kdlist = logspace(1,5,21);
kdlist = kdlist(1:9);
save("kd_rtot",'rel_eff_mat','rtotlist','kdlist')

rel_eff_mat = NaN(42,1);
for ii = 1:length(rel_eff_mat)
    fname = strcat("part2/tissue_300by900_opt",num2str(ii),".mat");
    if exist(fname,'file')~=0
        load(fname,"optMI","unifMI")
        rel_eff_mat(ii) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    end
end
rtotlist = logspace(2,5,7);
kdlist = logspace(-2,0.5,6);
save("kd_rtot_2",'rel_eff_mat','rtotlist','kdlist')

fnamelist = ["kd_rtot_2","kd_rtot"];
make_panel_7B(fnamelist);
% b = bar([unifreal,locreal]);
% b.FaceColor = 'flat';
% b.CData(6:end,:) = repmat([1 0 0],5,1);



