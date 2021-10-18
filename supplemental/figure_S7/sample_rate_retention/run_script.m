clear
%% setting paramaters
load('tissue_300by900_szopt','receptor_params')
load('scheme_parameter')

% setting feedback scheme parameter
thour = 2;
scheme_param.mean_cell_radius = ellipsePerimeter(2,5)/(2*pi);
conversion_factor = conc2count(scheme_param.mean_cell_radius,scheme_param.N);
scheme_param.T = 60*60*thour; % seconds
scheme_param.kd = receptor_params.kd*conversion_factor;
scheme_param.rtot = receptor_params.rtot;
scheme_param.receptornoise = receptor_params.receptornoise;
scheme_param.h = scheme_param.h*1.5;

% sampling rates to loop over
nparam = 7;
samplist = round(logspace(0,3,nparam));

schemelist = zeros(nparam,3);
uniflist = zeros(nparam,3);
statlist = zeros(nparam,4,3);
fnamelist = ["tissue_retention_1","tissue_retention_2","tissue_retention_3"];
for jj = 1:length(fnamelist)
    fname = fnamelist(jj);
    for ii = 1:nparam
        scheme_param.nsamp = samplist(ii);
        [schemerate,unifrate,statsummary] = racing_cells(fname,scheme_param,...
            "task","retention");
        files = dir(strcat(fname,"*"));
        dd = zeros(length(files));
        for j = 1:length(dd)
            dd(j) =datenum(files(j).date);
        end
        [~, id]=max(dd);
        curname = files(id).name;
        rename = strcat(curname(1:end-4),"_",num2str(ii),".mat");
        movefile(curname, rename);
%         load(strcat(fname,"_retention_feedback_",num2str(ii)));
        statlist(ii,:,jj) = statsummary;
        schemelist(ii,jj) = schemerate;
        uniflist(ii,jj) = unifrate;
    end
end
save("nsamp_retention_scan")

%% plotting
load("nsamp_retention_scan")
loglog(samplist, mean(schemelist,2),'linewidth',2,'Color',[0.8500 0.3250 0.0980]*1.12)
hold on
loglog(samplist, mean(uniflist,2),'linewidth',2,'Color',[0 0.4470 0.7410]*1.12)
hold off
xlabel('Sampling rate (per second)')
ylabel('Error rate (%)')
set(gca,'fontsize',27)
legend("Feedback scheme","Uniform",...
    'location','southeast','fontsize',25)
pbaspect([1,1,1])
legend("boxoff")
saveas(gca, "nsamp_scan_retention_plot.svg")