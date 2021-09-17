function [] = make_panel_5D_5E(fnamelist,fnamelist_grad)
%% PANEL D-E 
%data panel D (tissue)
telapsed = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'unif_time','scheme_time');
    telapsed = [telapsed;[unif_time,scheme_time]];
end


%data panel E (gradient)
telapsed_grad = [];
for ii = 1:length(fnamelist_grad)
    load(fnamelist_grad(ii),'unif_time','scheme_time');
    telapsed_grad = [telapsed_grad;[unif_time,scheme_time]];
end

% plots
palette = ["#297ec2","#f15e22"];
tiledlayout(2,2,'tilespacing','compact')
nexttile
for ii = 1:2
histogram(telapsed(:,ii),'Normalization','probability','BinWidth',8,...
    'FaceColor',palette(ii),'FaceAlpha',1,'linewidth',1)
hold on
end
hold off
xlim([0,6*64])
box off
set(gca,'fontsize',25)
pbaspect([2,1,1])

nexttile
disp(sum(telapsed<60)/size(telapsed,1)*100);
b = bar(sum(telapsed<60)/size(telapsed,1)*100,'linewidth',1);
b.FaceColor = 'flat';
b.CData(1,:) = [41, 126, 194]/255;
b.CData(2,:) = [241, 94, 34]/255;
ylim([0,100])
set(gca,'fontsize',25,'xticklabel',[])
box off
pbaspect([0.5,2,1])

nexttile %panel E histrogram
for ii = 1:2
histogram(telapsed_grad(:,ii),'Normalization','probability',...
        'BinWidth',8,'FaceColor',palette(ii),'FaceAlpha',1,'linewidth',1)
hold on
end
hold off
xlim([0,6*64])
box off
set(gca,'fontsize',25)
pbaspect([2,1,1])

nexttile
disp(sum(telapsed_grad<60)/size(telapsed_grad,1)*100);
b = bar(sum(telapsed_grad<60)/size(telapsed_grad,1)*100,'linewidth',1);
b.FaceColor = 'flat';
b.CData(1,:) = [41, 126, 194]/255;
b.CData(2,:) = [241, 94, 34]/255;
ylim([0,100])
set(gca,'fontsize',25,'xticklabel',[])
box off
pbaspect([0.5,2,1])
saveas(gcf,"panel_5D_5E.svg")

end