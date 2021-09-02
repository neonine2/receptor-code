function [ecmpts] = generateECM(N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% N = 1200*12*6;
% sz = gamrnd(75/3.7,3.7,1,N); %in um
sz = 75+5*randn(1,N);
X1 = 2400*rand(1,N);
Y1 = 1200*rand(1,N);
angle = 2*pi.*rand(1,N);
X2 = X1 + sz.*cos(angle);
Y2 = Y1 + sz.*sin(angle);
plot([X1; X2], [Y1; Y2],'k')
ecmpts = zeros(2000000,2);
lenspt=1;
for ii = 1:N
    segment = 0:0.5:sz(ii);
    pts = [X1(ii)+segment.*cos(angle(ii));Y1(ii)+segment.*sin(angle(ii))]';
    ind2 = lenspt+size(pts,1)-1;
    ecmpts(lenspt:ind2,:) = pts;
    lenspt = lenspt + size(pts,1);
end

end

