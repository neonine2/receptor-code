function [score] = clusterScore(rvec,c)
%CLUSTERSCORE Summary of this function goes here
nhbd = 10;

[~,I] = max(c);
m = length(rvec);
rvec = circshift(rvec,floor(m/2)-I);
c = circshift(c,floor(m/2)-I);

[~,I] = max(c);
score = sum(rvec((I-nhbd):(I+nhbd)))/((2*10+1)/m);
 
end

