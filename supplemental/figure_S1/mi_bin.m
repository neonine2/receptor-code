len = 50;
nbin = logspace(1,3,len);
nreceptor = 5;
receptor = logspace(0,2,nreceptor);
milist = zeros(len,nreceptor);
c=1;
for jj = 1:nreceptor
    for ii = 1:len
        m = nbin(ii);
        rvec = 1/m;
        params = struct('rtot',receptor(jj)*m,...
                            'kd',10,'receptornoise',0);
        mu = c/m;
        milist(ii,jj) = single_MI(rvec,mu,params)*m*log2(exp(1));
    end
end

% tiledlayout(2,1)
% nexttile
% semilogx(nbin, milist,'linewidth',1)
% set(gca,'fontsize',13,'linewidth',1)
% ylabel('MI')
% xlabel('number of bins')
% nexttile
n=256;
colors = [0.008,0.247,0.647;
          0.416,0.463,0.698;
          0.631,0.651,0.784;
          0.796,0.804,0.851;
          0.886,0.886,0.886];

hold on
for k=1:nreceptor
    semilogx(nbin, milist(:,k),'linewidth',1,'Color',colors(6-k,:),...
                                            'linewidth',2);
end
semilogx(nbin,c*(1-log(c)+log(nbin))*log2(exp(1)),'r','linewidth',2)
hold off
ylabel('I(C;A) (bits)')
xlabel('number of bins (m)')
legend('N/m = 1','N/m = 3.2','N/m = 10',...
       'N/m = 32','N/m = 100','Upper bound',...
                            'Location','northwest','box','off')
set(gca,'fontsize',17,'linewidth',1,"Xscale","log")
pbaspect([1,1,1])
