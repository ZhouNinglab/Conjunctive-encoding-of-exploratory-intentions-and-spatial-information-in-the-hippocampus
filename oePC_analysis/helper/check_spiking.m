function [cellx_spchecked] = check_spiking(binActivity,CellDif,explap,ObjBin)
%CHECK_SPIKING 此处显示有关此函数的摘要
%   此处显示详细说明

neuronBinLapS = binActivity.neuronBinLapS;
neuronBinL = CellDif.neuronBinL;
cellX1 = CellDif.cellXid;
n1 = numel(explap.explap1);
cellx_spchecked = [];
for j =  1:numel(cellX1)
    act1 = neuronBinLapS(cellX1(j),ObjBin,explap.explap1);
    act1 = transpose(reshape(act1,[numel(ObjBin),n1]));
    dif1 = neuronBinL(cellX1(j),:);
    
    sumtemp = sum(act1>0);
    if sum(sumtemp(find(dif1==1)) >=2) > 0
        cellx_spchecked = [cellx_spchecked cellX1(j)];
    end
   
end




end

