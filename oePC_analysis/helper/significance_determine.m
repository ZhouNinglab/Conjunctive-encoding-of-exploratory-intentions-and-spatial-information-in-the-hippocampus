function [Cells] = significance_determine(DataOriginal,DataShf,thr)
%SIGNIFICANCE_DETERMINE 此处显示有关此函数的摘要
%   DataOriginal mustbe n * 1 matrix
%   DataShf mustbe n * repeatNum matrix
%   Default threshold 1.65

    if nargin <= 2
        thr = 1.65;
    end
    d1 =  ndims(DataShf);
    data_mean = mean(DataShf,d1);
    data_std = std(DataShf,[],d1);
    
    
    SigVal = (DataOriginal - data_mean) ./ data_std;
    Cells.ZValue = SigVal;
    Cells.Determinant = (SigVal >= thr);
    Cells.ID = find(Cells.Determinant==1);
    Cells.ratio = numel(Cells.ID) / size(data_mean,1);
end

