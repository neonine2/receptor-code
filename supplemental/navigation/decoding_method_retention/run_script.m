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
scheme_param.nsamp = 30;

fnamelist = ["tissue_retention_1","tissue_retention_2","tissue_retention_3"];
nenv = length(fnamelist);
decodemethod = ["randomwalk","perfect","memory"];
nmethod = length(decodemethod);

schemelist = zeros(nmethod,nenv);
uniflist = zeros(nmethod,nenv);
statlist = zeros(nmethod,4,nenv);
for jj = 1:nenv
    fname = fnamelist(jj);
    for ii = 1:nmethod
        if isequal(decodemethod(ii), "memory")
            [schemerate,unifrate,statsummary] = ...
            racing_cells(fname,scheme_param,...
            "hasMemory",true,"task","retention");
        else
            [schemerate,unifrate,statsummary] = ...
            racing_cells(fname,scheme_param,...
            "decoder_method",decodemethod(ii),"task","retention");
        end
        statlist(ii,:,jj) = statsummary;
        schemelist(ii,jj) = schemerate;
        uniflist(ii,jj) = unifrate;
    end
end
save("decodemethod_retention_scan")

%% plotting
clear
fnamelist = ["tissue_retention_1_retention_feedback_randomwalk",...
              "tissue_retention_2_retention_feedback_randomwalk",...
                "tissue_retention_3_retention_feedback_randomwalk"];
telapsed_unif = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'unif_traj','scheme_traj');
    telapsed_unif = [telapsed_unif;[unif_traj(:),scheme_traj(:)]+1];
end
fnamelist = ["tissue_retention_1_retention_feedback_perfect",...
              "tissue_retention_2_retention_feedback_perfect",...
                "tissue_retention_3_retention_feedback_perfect"];
telapsed = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'unif_traj','scheme_traj');
    telapsed = [telapsed;[unif_traj(:),scheme_traj(:)]+1];
end
fnamelist = ["tissue_retention_1_retention_feedback_Memory",...
              "tissue_retention_2_retention_feedback_Memory",...
                "tissue_retention_3_retention_feedback_Memory"];
telapsed_memory = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'unif_traj','scheme_traj');
    telapsed_memory = [telapsed_memory;[unif_traj(:),scheme_traj(:)]+1];
end

tiledlayout(1,4,'tilespacing','compact')
% plots
palette = ["#297ec2","#f15e22"];
nexttile
for ii = 1:2
histogram(telapsed_unif(:,ii),'Normalization','probability',...
    'BinWidth',1,...
    'FaceColor',palette(ii),'FaceAlpha',1)
hold on
end
hold off
xlim([0,30])
ylim([0,1])
box off
set(gca,'fontsize',12)
ylabel('Proportion of time steps')
xlabel('Distance from source')
pbaspect([1,1,1])

nexttile
for ii = 1:2
histogram(telapsed(:,ii),'Normalization','probability',...
    'BinWidth',1,...
    'FaceColor',palette(ii),'FaceAlpha',1)
hold on
end
hold off
xlim([0,15])
box off
set(gca,'fontsize',12)
ylabel('Proportion of time steps')
xlabel('Distance from source')
pbaspect([1,1,1])

nexttile
for ii = 1:2
histogram(telapsed_memory(:,ii),'Normalization','probability',...
    'BinWidth',1,...
    'FaceColor',palette(ii),'FaceAlpha',1)
hold on
end
hold off
xlim([0,15])
box off
set(gca,'fontsize',12)
ylabel('Proportion of time steps')
xlabel('Distance from source')
pbaspect([1,1,1])

load("decodemethod_retention_scan")
nexttile
X = categorical({'Random','Maximal increase','Temporal averaging'});
X = reordercats(X,{'Random','Maximal increase','Temporal averaging'});
b = bar(X,[mean(uniflist');mean(schemelist')]');
b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(1).CData = repmat([41, 126, 194]/255,3,1);
b(2).CData = repmat([241, 94, 34]/255,3,1);
ylim([0,100])
ylabel('Error rate (%)')
set(gca,'fontsize',12)
box off
pbaspect([1,1,1])
saveas(gca,"decodermethod_retention.svg")

