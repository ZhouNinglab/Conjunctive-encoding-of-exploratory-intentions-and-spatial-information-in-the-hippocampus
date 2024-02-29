function [bsResult] = bootstrap_analysis(data,n1,repeatNum)
%BOOTSTRAP_ANALYSIS 此处显示有关此函数的摘要
%   此处显示详细说明
% if ndims(data) > 1
%     disp('Input data for bootstrap analysis must have only 1 dimension!')
%     return
% end

if nargin < 3
    repeatNum = 1000;    
end


bsResult.repeatNum = repeatNum;
mean_ori = group_mean(data,n1);
bsResult.mean_ori = mean_ori;

 mean_b = bootstrp(repeatNum,@(data)group_mean(data,n1),data);


bsResult.mean_b = mean_b;
ratio = sum(mean_b < mean_ori) / repeatNum;
bsResult.ratio = ratio;

end

