function [meanDif ] = group_mean(data,n1)
%GROUP_MEAN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

mean1 = mean(data(1:n1));
mean2 = mean(data(n1+1:end));
meanDif = mean1 - mean2;



end

