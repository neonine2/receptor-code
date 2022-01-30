clear
mode = "diff";
rtotlist = logspace(2,5,13);
kdlist = logspace(-2,2.5,19);
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
        "index",jj,'optgrid',[5,15],"mode",mode);
end 
save("kd_rtot_concfixed");

% plotting
mode = "diff";
rel_eff_mat = NaN(13*19,1);
for ii = 1:length(rel_eff_mat)
    fname = strcat("tissue_300by900_opt",num2str(ii),".mat");
    if exist(fname,'file')~=0
        load(fname,"optMI","unifMI")
        rel_eff_mat(ii) = compute_efficacy(mean(optMI),mean(unifMI),"mode",mode);
    end
end
rtotlist = logspace(2,5,13);
kdlist = logspace(-2,2.5,19);
save("kd_rtot",'rel_eff_mat','rtotlist','kdlist')
fnamelist = ["kd_rtot"];
[locreal, unifreal] = make_panel_7B(fnamelist);
b = bar([unifreal,locreal]);
b.FaceColor = 'flat';
b.CData(6:end,:) = repmat([1 0 0],5,1);

