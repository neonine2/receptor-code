clear
%% setting paramaters
load('tissue_300by900_szopt','receptor_params')
load('scheme_parameter')

% setting feedback scheme parameter
thour = 1;
conversion_factor = conc2count(scheme_param.mean_cell_radius,scheme_param.N);
scheme_param.T = 60*60*thour; % seconds
scheme_param.kd = receptor_params.kd*conversion_factor;
scheme_param.rtot = receptor_params.rtot;
scheme_param.receptornoise = receptor_params.receptornoise;

% default parameter
fnamelist = ["tissue_300by900",...
             "tissue_taxis_1",...
             "tissue_taxis_2"];
gammalist = [0,2e-2,1e-1];
         
for jj = [2,3]
    scheme_param.gamma = gammalist(jj); % for w1dist
    for ii = 2:3
        fname = fnamelist(ii);
        racing_cells(fname, scheme_param,...
            "envmodel","tissue",...
            "receptor","w1dist",...
            "task","localization",...
            "verbose",true);
    end
end

%%plotting
clear
filelist = dir("tissue*");
for ii = 1:3
    load(filelist(ii).name,"schemerate","unifrate")
    schemedata(ii) = schemerate;
    unifdata(ii) = unifrate;
end
[schemedata(2:3)] = [schemedata(3),schemedata(2)];
unifdata = mean(unifdata)*ones(1,3);

b = bar([unifdata;schemedata]');
b(1).CData = [41, 126, 194]/255;
b(2).CData = [241, 94, 34]/255;
set(gca,'fontsize',14,'xticklabel',["0","0.02","0.1"])
xlabel("$\gamma$",'interpreter','latex','fontsize',22)
ylabel("Success rate (%)",'fontsize',16)
box off
legend("Uniform","Dynamic protocol",'Location','northwest')
legend("boxoff")
pbaspect([2,1,1])
saveas(gca,"w1dist_localization.svg")