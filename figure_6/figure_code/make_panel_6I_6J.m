function [] = make_panel_6I_6J(fnamelist,fnamelist_grad)

%% PANEL I-J
telapsed = [];
for ii = 1:length(fnamelist)
    load(fnamelist(ii),'unif_traj','scheme_traj');
    telapsed = [telapsed;[unif_traj(:),scheme_traj(:)]-1+2];
end


telapsed_grad = [];
for ii = 1:length(fnamelist_grad)
    load(fnamelist_grad(ii),'unif_traj','scheme_traj');
    telapsed_grad = [telapsed_grad;[unif_traj(:),scheme_traj(:)]-1+2];
end

% plots
palette = ["#297ec2","#f15e22"];
tiledlayout(2,2)
nexttile
for ii = 1:2
histogram(telapsed(:,ii),'Normalization','probability',...
    'BinWidth',0.2,'FaceColor',palette(ii),'FaceAlpha',1)
hold on
end
hold off
xlim([0,10])
box off
set(gca,'fontsize',25)
pbaspect([2.5,1,1])

nexttile
disp(strcat(num2str(sum(telapsed>5)),"/",...
    num2str(size(telapsed,1))))
b = bar(sum(telapsed>5)/size(telapsed,1)*100);
b.FaceColor = 'flat';
b.CData(1,:) = [41, 126, 194]/255;
b.CData(2,:) = [241, 94, 34]/255;
ylim([0,17])
set(gca,'fontsize',25,'xticklabel',[])
box off
pbaspect([0.5,2,1])

nexttile %panel J histrogram
for ii = 1:2
histogram(telapsed_grad(:,ii),'Normalization','probability',...
        'BinWidth',0.2,'FaceColor',palette(ii),'FaceAlpha',1)
hold on
end
hold off
xlim([0,10])
box off
set(gca,'fontsize',25)
pbaspect([2.5,1,1])

nexttile
disp(strcat(num2str(sum(telapsed_grad>5)),"/",...
    num2str(size(telapsed_grad,1))))
b = bar(sum(telapsed_grad>5)/size(telapsed_grad,1)*100);
b.FaceColor = 'flat';
b.CData(1,:) = [41, 126, 194]/255;
b.CData(2,:) = [241, 94, 34]/255;
ylim([0,17])
set(gca,'fontsize',25,'xticklabel',[])
box off
pbaspect([0.5,2,1])
saveas(gcf,"panel_6I_6J.svg")

end