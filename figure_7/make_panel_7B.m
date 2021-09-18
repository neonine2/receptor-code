function [] = make_panel_7B(fnamelist)
load(fnamelist(2),'rel_eff_mat','rtotlist','kdlist')
rel_eff = rel_eff_mat;
paramlist = combvec(rtotlist,kdlist)';

%interpolation
F = scatteredInterpolant(paramlist(:,1),paramlist(:,2),...
                                            rel_eff);
rtotlist = logspace(2,5,16);
kdlist = logspace(-2,0.5,12);
paramlist = combvec(rtotlist,kdlist)';
rel_eff_FULL = F(paramlist(:,1),paramlist(:,2));

load(fnamelist(1),'rel_eff_mat','rtotlist','kdlist')
rel_eff = cat(1,rel_eff_FULL, rel_eff_mat);
paramlist = cat(1,paramlist, combvec(rtotlist,kdlist)');
F = scatteredInterpolant(paramlist(:,1),paramlist(:,2),...
                                            rel_eff);
rtotlist = logspace(2,5,42*4);
kdlist = logspace(-2,5,62*4);
paramlist = combvec(rtotlist,kdlist)';
rel_eff_FULL = reshape(F(paramlist(:,1),paramlist(:,2)),...
    length(rtotlist), length(kdlist))';

c = viridis;
colormap(c);
% rel_eff_FULL((rel_eff_FULL<10.3) .* (rel_eff_FULL>9.8)==1) = 200;
imagesc(rtotlist,kdlist,rel_eff_FULL) % full heatmap
set(gca,'XScale','log','YScale','log','fontsize',26);
pbaspect([1,2,1])
colorbar()

%% adding empirical value of real receptors
hold on
scatter([32000,1572,593,200,22300],...
            [25000,5,4,12,235],70,'red','filled',...
                            'MarkerEdgeColor','k')
scatter([1000,1000,10000,17000,50000],...
            [0.013,0.019,0.05,12.5,2],70,'white',...
                            'filled','MarkerEdgeColor','k')
hold off
saveas(gca,"panel_7B.svg")
end