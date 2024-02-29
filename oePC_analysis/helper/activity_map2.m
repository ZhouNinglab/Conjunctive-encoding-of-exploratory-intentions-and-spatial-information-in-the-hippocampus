function [activitymap2] = activity_map2(activitymap1,opt)
%ACTIVITY_MAP2 此处显示有关此函数的摘要
%   Generate mean activity maps for multiple neurons, with peaks in
%   temporal sequence
%   Opt == 0, peak values are normalized to 1
%   Opt == 1, peak values are not normalized


[M1,B1] = max(activitymap1,[],2);

if nargin < 2
    opt = 1;
end

if  opt == 0;
    map2 = activitymap1;
elseif  opt == 1;
    
    map2 = activitymap1./M1;
end

B2 = [B1 (1:numel(B1))'];
B2 = sortrows(B2, 1, 'ascend');
idx = B2(:,2);
activitymap2.map = map2(idx,:);
activitymap2.idx = idx;


end

