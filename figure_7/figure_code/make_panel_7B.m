function [loceff,unifeff,eff_FULL] = make_panel_7B(fnamelist)

ub=-1;lb=-1;
param = [];
eff = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'rel_eff_mat','rtotlist','kdlist')
    paramlist = combvec(rtotlist,kdlist)';
    %get only param with valid efficacy
    paramlist = paramlist(~isnan(rel_eff_mat),:);
    rel_eff_mat = rel_eff_mat(~isnan(rel_eff_mat));
    outbnd = prod(paramlist<=ub & paramlist>=lb,2);
    paramlist = paramlist(~outbnd,:);
    rel_eff_mat = rel_eff_mat(~outbnd);
    ub = log10(max(paramlist));
    lb = log10(min(paramlist));
    
    %fill in blanks with interpolant
    F = scatteredInterpolant(paramlist(:,1),paramlist(:,2),...
                                        rel_eff_mat);
    rtotlist = logspace(lb(1),ub(1),(ub(1)-lb(1))*20);
    kdlist = logspace(lb(2),ub(2),(ub(2)-lb(2))*20);
    paramlist = combvec(rtotlist,kdlist)';
    rel_eff_mat = F(paramlist(:,1),paramlist(:,2));

    param = [param;paramlist];
    eff = [eff;rel_eff_mat];
        
    %set new limit to remove overlapping points from other files
    ub = max(param);
    lb = min(param);
end
param = param(eff>=0,:);
eff = eff(eff>=0);
ub = log10(max(param));
lb = log10(min(param));
rtotlist = logspace(lb(1),ub(1),(ub(1)-lb(1))*60);
kdlist = logspace(lb(2),ub(2),(ub(2)-lb(2))*60);
paramlist = combvec(rtotlist,kdlist)';

F = scatteredInterpolant(param(:,1),param(:,2),eff);
eff_FULL = reshape(F(paramlist(:,1),paramlist(:,2)),...
                    length(rtotlist), length(kdlist))';

% Define rgb matrix first (matrix size ncolors x 3)
% cmap = [0.945,0.945,0.945;
%           0.784,0.886,0.922;
%           0.655,0.776,0.867;
%           0.553,0.639,0.792;
%           0.486,0.482,0.698;
%           0.447,0.306,0.584;
%           0.420,0.000,0.467];
% cmap = interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,size(cmap,1)*4));
cmap = "Turbo";

% colormap(cmap);
% rel_eff_FULL((rel_eff_FULL<10.3) .* (rel_eff_FULL>9.8)==1) = 200;

colormap(cmap);
imagesc(rtotlist,kdlist,eff_FULL) % full heatmap
set(gca,'XScale','log','YScale','log','fontsize',20);
pbaspect([2,4,1])
%     cmocean('-tempo');
colorbar()

%% adding empirical value of real receptors
hold on
% scatter([32000,1572,593,200,22300],...
%             [25000,5,4,12,235],70,'red','filled',...
%                             'MarkerEdgeColor','k')
locN = [1100,1572,593,200,22300];
lockd = [1.53,5,4,12,235];
scatter(locN,lockd,70,'red','filled','MarkerEdgeColor','k')
unifN = [1000,1000,10000,40000,50000];
unifkd = [0.013,0.019,0.05,12,2];
scatter(unifN,unifkd,70,'white','filled','MarkerEdgeColor','k')
hold off
loceff = F(locN,lockd);
unifeff = F(unifN,unifkd);
saveas(gca,"panel_7B.svg")
end