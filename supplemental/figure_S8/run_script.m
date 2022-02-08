clear all;
close all;

% testing efficiency of feedback scheme in both static and dynamic
% environment and across different parameters
fname = ["tissue_300by900_szopt.mat","soil_var_2_szopt.mat"]; %tissue, soil
cellmodel = 2; %10um radius cell

% load parameters for feedback scheme simulation
load('scheme_parameter','scheme_param') 
ind = [19,20];

for ii = 1:length(fname)
    feedback_scheme_efficiency(which(fname(ii)),scheme_param,...
                          "dual","static",cellmodel,"selectind",ind(ii));
end

%plotting
fname = ["tissue_300by900_szopt_scheme_dual_static_10um",...
             "soil_var_2_szopt_scheme_dual_static_10um"];
hlist = logspace(-3,-1,16);
kofflist = logspace(-2,0,16);

tiledlayout(2,1)
for ii = 1:2
    load(fname(ii),'schemer','schemec','paramlist','recstat');
    nparam = size(schemer,3);
    score = zeros(nparam,1);
    for jj = 1:nparam
        score(jj) = clusterScore(schemer(:,:,jj),schemec);
    end
    score = score(~isnan(score));
    h = squeeze(recstat(:,1,~isnan(score)));
    koff = paramlist(~isnan(score),2);
    F = scatteredInterpolant(h,koff,score);
    
    newha = logspace(log10(min(recstat(:,1,:))+2e-5),...
                   log10(max(recstat(:,1,:))-1e-5),length(hlist));
    newkoff = logspace(log10(min(kofflist)),...
                   log10(max(kofflist)),length(kofflist));
    newcoord = combvec(newha,newkoff)';
    score = F(newcoord(:,1),newcoord(:,2));
    score = reshape(score,length(newha),length(newkoff));
    nexttile
    colormap(viridis)
    imagesc(newkoff, newha, score)
    set(gca,'XScale','log','YScale','log','fontsize',13);
    caxis([0.006,4.87])
    caxis
    xlim([0.01,1])
    colorbar()
    pbaspect([1,1,1])
    xlabel('$k_{\mathrm{off}}$','interpreter','latex')
    ylabel('$\langle h A_i \rangle$','interpreter','latex')
end

%% plotting call out boxes 
clear all;
close all;

fname = ["tissue_300by900_szopt_scheme_dual_static_10um",...
             "soil_var_2_szopt_scheme_dual_static_10um"];
[envtissue, envsoil] = make_panel_5B(["tissue_300by900_szopt_scheme_koff_static_10um",...
                     "soil_var_2_szopt_scheme_koff_static_10um"]);
                 
% Tissue
paramindex = [65,98,120,160-1];
tiledlayout(2,2,'TileSpacing','compact')
for ii = 1:4
    load(fname(1),'schemer','scheme_param','recstat','paramlist');
    nexttile
    colormap(gray)
    imagesc(envtissue)
    pbaspect([1,1,1])
    m = scheme_param.N;
    r = schemer(1,:,paramindex(ii))';
    angle = linspace(0,2*pi*(1-1/m),m);
    s = scheme_param.mean_cell_radius*10;
    cellsurf = 221 + (s+r*m*5).*[cos(angle)',sin(angle)'];
    hold on
    fill(cellsurf(:,1),cellsurf(:,2),[17 17 17]/18,'Linewidth',1)
    rectangle('Position',[221-s,221-s,2*s,2*s],'Curvature',[1,1], ...
        'FaceColor', [17 17 17]/32, 'EdgeColor', 'k','Linewidth',2);
    title(recstat(1,1:2,paramindex(ii)))
    set(gca,'xticklabels',[],'yticklabels',[],'xtick',[],'ytick',[],...
        'linewidth',2);
%     colorbar('Ticks',[0.5,1,1.5],'TickLabels',["0.5","1","1.5"],...
%         'fontsize',12)
end

% Soil
paramindex = [65,98,120,160-1];
tiledlayout(2,2,'TileSpacing','compact')
for ii = 1:4
    load(fname(2),'schemer','scheme_param','recstat','paramlist');
    nexttile
    colormap(gray)
    imagesc(envsoil)
    pbaspect([1,1,1])
    m = scheme_param.N;
    r = schemer(1,:,paramindex(ii))';
    angle = linspace(0,2*pi*(1-1/m),m);
    s = scheme_param.mean_cell_radius*10;
    cellsurf = 221 + (s+r*m*5).*[cos(angle)',sin(angle)'];
    hold on
    fill(cellsurf(:,1),cellsurf(:,2),[17 17 17]/18,'Linewidth',1)
    rectangle('Position',[221-s,221-s,2*s,2*s],'Curvature',[1,1], ...
        'FaceColor', [17 17 17]/32, 'EdgeColor', 'k','Linewidth',2);
    title(recstat(1,1:2,paramindex(ii)))
    set(gca,'xticklabels',[],'yticklabels',[],'xtick',[],'ytick',[],...
        'linewidth',2);
%     colorbar('Ticks',[0.5,1,1.5],'TickLabels',["0.5","1","1.5"],...
%         'fontsize',12)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

fname = ["tissue_300by900_szopt.mat","soil_var_2_szopt.mat"];
testing_param = ["koff","h"];
mode = ["static","dynamic"];
cellmodel = 1; %5um
% load parameters for feedback scheme simulation
load('scheme_parameter','scheme_param') 

for ii = 1:length(fname)
    for jj = 1:length(testing_param)
        for kk = 1:length(mode)
            feedback_scheme_efficiency(which(fname(ii)),scheme_param,...
                testing_param(jj),mode(kk),cellmodel)
        end
    end
end

%% plotting

% IMPORTANT: make sure to rename file as this will generate panel5 name
% file
make_panel_5C(["tissue_300by900_szopt_scheme_h_static_5um",...
                "soil_var_2_szopt_scheme_h_static_5um",...
                  "tissue_300by900_szopt_scheme_koff_static_5um",...
                     "soil_var_2_szopt_scheme_koff_static_5um"])
                 
make_panel_5F(["tissue_300by900_szopt_scheme_h_dynamic_5um",...
                     "soil_var_2_szopt_scheme_h_dynamic_5um",...
                     "tissue_300by900_szopt_scheme_koff_dynamic_5um",...
                     "soil_var_2_szopt_scheme_koff_dynamic_5um"])