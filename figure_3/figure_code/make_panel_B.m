function [] = make_panel_B(optfname, perturbfname, varargin)

p = inputParser;
addRequired(p,'optfname',@isfile);
addRequired(p,'perturbfname',@isfile);

% optional arguments
addParameter(p,'save_plot',true,@islogical);
parse(p,optfname,perturbfname,varargin{:});
optfname = p.Results.optfname;
perturbfname = p.Results.perturbfname;
save_plot = p.Results.save_plot;

%% loading data
load(perturbfname)
load(optfname,'unifMI', 'optMI', 'optr', 'envmean')
[nloc,m] = size(optr,[1,2]);
unifMI = unifMI(:,cellmodel);
optr = optr(:,:,cellmodel);
optMI = optMI(:,cellmodel);
envmean = envmean(:,:,cellmodel);
opt_eff = (sum(optMI) - sum(unifMI))./sum(unifMI) .* 100;

%% panel B heatmap
perturbd_eff = zeros(length(shiftparam)*length(flattenparam),1);
for ii = 1:nparam
    pMI = perturbMI(:,ii);
    perturbd_eff(ii) = (sum(pMI) - sum(unifMI))./sum(unifMI) .* 100;
end
F = scatteredInterpolant(perturbparam,perturbd_eff,'natural','linear');
shiftp = linspace(-15,15,81);
flattenp = linspace(1,25,81);
newparam = combvec(shiftp,flattenp)';
interpmat = reshape(F(newparam),length(shiftp),length(flattenp));
interpmat = interpmat./opt_eff;
disp(['Max = ', num2str(max(perturbd_eff,[],'all')), ' Opt = ', num2str(opt_eff)])

%plotting
figure(1)
c = viridis;
colormap(c(35:end,:));
imagesc(flattenp, shiftp/m*360, interpmat)
yticks([-45,-30,-15,0,15,30,45])
set(gca,'fontsize',20)
colorbar('location','eastoutside','ticks',0:0.2:1, 'LineWidth',0.5)
pbaspect([1 1 1])

if save_plot
    [filepath,name,ext] = fileparts(perturbfname);
    saveas(gcf,strcat(name,"_heatmap"),'svg')
end


%% panel B receptor callout
figure(2)
sublist = [11,4;9,6;3,7];
tiledlayout(1,size(sublist,1),'tilespacing','compact')
for ii = 1:size(sublist,1)
    param_sub = sublist(ii,:);
    shifted_perturbr = zeros(nloc,m);
    shifted_optr = zeros(nloc,m);
    ind = sub2ind([length(shiftparam),length(flattenparam)],...
        param_sub(1),param_sub(2));
    disp(['xtrans = ',num2str(perturbparam(ind,1)/m*360),...
        ' flatf = ',num2str(perturbparam(ind,2)),...
        ' %DeltaI/I = ',num2str(perturbd_eff(ind)/opt_eff*100)]);
    for jj = 1:nloc
        [~,I] = max(envmean(jj,:));
        shifted_perturbr(jj,:) = circshift(perturbr(jj,:,ind),m/2-I);
        shifted_optr(jj,:) = circshift(optr(jj,:),m/2-I);
    end
    
    mean_optr = mean(shifted_optr); % optimal mean r
    std_dev = std(shifted_optr);
    range_optr = [mean_optr + std_dev, fliplr(mean_optr - std_dev)];
    mean_perturbr = mean(shifted_perturbr); % shifted mean r
    std_dev = std(shifted_perturbr);
    range_perturbr = [mean_perturbr + std_dev, fliplr(mean_perturbr - std_dev)];
    
    nexttile
    fill([1:m, fliplr(1:m)],range_optr,[1,1,1]/1.5,'LineStyle','none');
    hold on
    plot(1:m, mean_optr,'color',[1,1,1]/2,'LineWidth',0.75);
    fill([1:m, fliplr(1:m)],range_perturbr,[1,0.8,0.8],'LineStyle','none');
    plot(1:m, mean_perturbr, 'r', 'LineWidth', 0.75);
    pbaspect([1,1.2,1])
    if ii == 1
        yticks([0.04,0.35])
        yticklabels({'4%','35%'})
    else
        set(gca,'YTicklabel',[]);
    end
    ylim([0,0.35])
    set(gca,'XTickLabel',[],'fontsize',18,'linewidth',0.3);
    hold off
end

if save_plot
    [filepath,name,ext] = fileparts(perturbfname);
    saveas(gcf,strcat(name,"_rvec"),'svg');
end
end