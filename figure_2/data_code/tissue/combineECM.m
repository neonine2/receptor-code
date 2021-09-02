m = 1;
ecmset = zeros(3000000,2,m);

for ii = 1:m
%     ecmset(:,:,ii) = generateECM(43200); %ecm_4, m=4
    ecm = generateECM(40000);
    ecmset(1:length(ecm),:,ii) = ecm;
    disp(ii);
end
ecmset = ecm;

ecmMAT = zeros(2200,1000);
for ii = 1:m
    [M,~,~] = histcounts2(ecmset(:,1,ii),ecmset(:,2,ii),100:2300,100:1100);
    ecmMAT = ecmMAT + M;
end

imagesc(ecmMAT')
colorbar()
% save('ecm_4.mat','ecmMAT','-v7.3')
save('ecm.mat','ecmMAT','-v7.3')
