function [meanDif ] = group_mean(data,n1)
%GROUP_MEAN 此处显示有关此函数的摘要
%   此处显示详细说明

mean1 = mean(data(1:n1));
mean2 = mean(data(n1+1:end));
meanDif = mean1 - mean2;



end

