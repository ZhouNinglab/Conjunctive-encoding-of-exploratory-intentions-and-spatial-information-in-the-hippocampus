function [ LapBin] = activity_map(binActivity,id,lapseq)
%ACTIVITY_MAP 此处显示有关此函数的摘要
%   Generate the activity map by laps for individual neurons

   neuronBinLapS = binActivity.neuronBinLapS;
   [n, nbin, lapnum] = size(neuronBinLapS);
   if nargin < 3
       lapseq = 1:lapnum; 
   end
   BinLap=neuronBinLapS(id,:,lapseq);
   LapBin=transpose(reshape(BinLap,nbin, []));
   

end

