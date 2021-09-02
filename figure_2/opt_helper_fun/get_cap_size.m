function [output, sub, ind] = get_cap_size(rvec, nhbd)
%get the proportion of receptors within +/- nhbd around peak
[M,I] = max(rvec);
m = length(rvec);
if sum(rvec==M) > 1
    warning('multiple maxima')
end

if nhbd*2+1 > m
    error('neighbourhood size too big')
end

if I < nhbd + 1
    ind = [1:I,(I+1):(I+nhbd),(m-(nhbd-I)):m];
elseif I > m - nhbd
    ind = [(I-nhbd):(I-1),I:m,1:(nhbd-(m-I))];
else
    ind = I-nhbd:I+nhbd;
end

sub = rvec(ind);
output = sum(sub);
 
end

