clear
%% setting paramaters
load('tissue_300by900_szopt','receptor_params')
load('scheme_parameter')

% setting feedback scheme parameter
thour = 6; % keep greater than or equal to 1
conversion_factor = conc2count(scheme_param.mean_cell_radius,scheme_param.N);
scheme_param.T = 60*60*thour; % seconds
scheme_param.kd = receptor_params.kd*conversion_factor;
scheme_param.rtot = receptor_params.rtot;
scheme_param.receptornoise = receptor_params.receptornoise;
scheme_param.nsamp = 30;

fnamelist = ["tissue_300by900","tissue_taxis_1",...
             "tissue_taxis_2"];
nenv = length(fnamelist);
decodemethod = ["randomwalk","perfect","memory"];
nmethod = length(decodemethod);

schemelist = zeros(nmethod,nenv);
uniflist = zeros(nmethod,nenv);
statlist = zeros(nmethod,4,nenv);
for jj = 1:nenv
    fname = fnamelist(jj);
    for ii = 3:nmethod
        if isequal(decodemethod(ii), "memory")
            [schemerate,unifrate,statsummary] = ...
            racing_cells(fname,scheme_param,"hasMemory",true);
        else
            [schemerate,unifrate,statsummary] = ...
            racing_cells(fname,scheme_param,"decoder_method",decodemethod(ii));
        end
        statlist(ii,:,jj) = statsummary;
        schemelist(ii,jj) = schemerate;
        uniflist(ii,jj) = unifrate;
    end
end
save("decodemethod_scan")

%% plotting
fnamelist = ["tissue_300by900_localization_feedback_randomwalk",...
              "tissue_taxis_1_localization_feedback_randomwalk",...
                "tissue_taxis_2_localization_feedback_randomwalk"];
telapsed_random = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'unif_time','scheme_time');
    telapsed_random = [telapsed_random;[unif_time,scheme_time]];
end
fnamelist = ["tissue_300by900_localization_feedback_perfect",...
              "tissue_taxis_1_localization_feedback_perfect",...
                "tissue_taxis_2_localization_feedback_perfect"];
telapsed = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'unif_time','scheme_time');
    telapsed = [telapsed;[unif_time,scheme_time]];
end
fnamelist = ["tissue_300by900_localization_feedback_Memory",...
              "tissue_taxis_1_localization_feedback_Memory",...
                "tissue_taxis_2_localization_feedback_Memory"];
telapsed_memory = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'unif_time','scheme_time');
    telapsed_memory = [telapsed_memory;[unif_time,scheme_time]];
end

tiledlayout(1,4,'tilespacing','compact')
% plots
palette = ["#297ec2","#f15e22"];
nexttile
for ii = 1:2
histogram(telapsed_random(:,ii),'Normalization','probability',...
    'BinWidth',8,...
    'FaceColor',palette(ii),'FaceAlpha',1)
hold on
end
hold off
xlim([0,6*64])
box off
set(gca,'fontsize',12)
ylabel('Proportion of cells')
xlabel('Time taken to source (min)')
pbaspect([1,1,1])

nexttile
for ii = 1:2
histogram(telapsed(:,ii),'Normalization','probability',...
    'BinWidth',8,...
    'FaceColor',palette(ii),'FaceAlpha',1)
hold on
end
hold off
xlim([0,6*64])
box off
set(gca,'fontsize',12)
ylabel('Proportion of cells')
xlabel('Time taken to source (min)')
pbaspect([1,1,1])

nexttile
for ii = 1:2
histogram(telapsed_memory(:,ii),'Normalization',...
    'probability','BinWidth',8,...
    'FaceColor',palette(ii),'FaceAlpha',1)
hold on
end
hold off
xlim([0,6*64])
box off
set(gca,'fontsize',12)
ylabel('Proportion of cells')
xlabel('Time taken to source (min)')
pbaspect([1,1,1])

load("decodemethod_scan")
nexttile
X = categorical({'Random','Maximal increase','Optimal with memory'});
X = reordercats(X,{'Random','Maximal increase','Optimal with memory'});
b = bar(X,[mean(uniflist');mean(schemelist')]');
b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(1).CData = repmat([41, 126, 194]/255,3,1);
b(2).CData = repmat([241, 94, 34]/255,3,1);
ylim([0,100])
set(gca,'fontsize',12)
ylabel('Success rate (%)')
box off
pbaspect([1,1,1])
saveas(gca,"decodermethod_localization.svg")

