clear
%% setting paramaters
load('tissue_300by900_szopt','receptor_params')
load('scheme_parameter')

% setting feedback scheme parameter
thour = 1; % keep greater than or equal to 1
conversion_factor = conc2count(scheme_param.mean_cell_radius,scheme_param.N);
scheme_param.T = 60*60*thour; % seconds
scheme_param.kd = receptor_params.kd*conversion_factor;
scheme_param.rtot = receptor_params.rtot;
scheme_param.receptornoise = receptor_params.receptornoise;
nparam = 7;
samplist = round(logspace(0,3,nparam));

schemelist = zeros(nparam,3);
uniflist = zeros(nparam,3);
statlist = zeros(nparam,4,3);
fnamelist = ["tissue_300by900","tissue_taxis_1",...
             "tissue_taxis_2"];
for jj = 1:length(fnamelist)
    fname = fnamelist(jj);
    for ii = 1:nparam
        scheme_param.nsamp = samplist(ii);
        [schemerate,unifrate,statsummary] = racing_cells(fname, scheme_param);
        files = dir(strcat(fname,"*"));
        dd = zeros(length(files));
        for j = 1:length(dd)
            dd(j) =datenum(files(j).date);
        end
        [~, id]=max(dd);
        curname = files(id).name;
        rename = strcat(curname(1:end-4),"_",num2str(ii),".mat");
        movefile(curname, rename);
        statlist(ii,:,jj) = statsummary;
        schemelist(ii,jj) = schemerate;
        uniflist(ii,jj) = unifrate;
    end
end
save("nsamp_scan")

%% plotting
load("nsamp_scan")
loglog(samplist, mean(schemelist,2),...
    'linewidth',1.5,'Color',[241, 94, 34]/255)
hold on
loglog(samplist, mean(uniflist,2),...
    'linewidth',1.5,'Color',[41, 126, 194]/255)
hold off
xlabel('Sampling rate (# samples per second)')
ylabel('Success rate (%)')
set(gca,'fontsize',16)
pbaspect([1.5,1,1])
legend("Feedback scheme","Uniform",...
    'location','southeast','fontsize',14)
legend("boxoff")
saveas(gca, "nsamp_scan_plot.svg")