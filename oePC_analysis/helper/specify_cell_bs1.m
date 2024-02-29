function [CellDif] = specify_cell(NeuronSL,explap,Lap)
%SPECIFY_LAP 此处显示有关此函数的摘要
%   此处显示详细说明

[n,m] = size(NeuronSL);
% neuronBinLapS = binActivity.neuronBinLapS;

% laps_theta_post = binActivity.laps_theta_post;
% [n,nbin, lapnum] = size(neuronBinLapS);

explap1 = explap.explap1;
explap2 = explap.explap2;


n1 = numel(explap1);
n2 = numel(explap2);
LapNew = [Lap(explap.explap1,:); Lap(explap.explap2,:)];
LapNum = size(LapNew,1);

neuronBinL = zeros(n,1);
neuronBinL = logical(neuronBinL);
data_all = zeros(n,LapNum);
cellXid = [];
fms = explap.fms;
fmsnew = fms([explap1; explap2]);
CellDif.Ratio = zeros(n,1);

for j1 = 1:n
    for j2 = 1:LapNum
 

%         data1 = neuronBinLapS(j1, ObjBin(j2),[explap1;explap2]);   
%         data1 = reshape(data1,[[n1+n2],1]);
        
%         data_all(j1,j2) = mean(double(NeuronSL(j1,fmsnew{j2})));
          data_all(j1,j2) = max(double(NeuronSL(j1,fmsnew{j2})),[],2);
   end
 
    data1 = data_all(j1,:);
    bsResult = bootstrap_analysis(data1,n1,10000);
    CellDif.Ratio(j1) = bsResult.ratio;
    if  (bsResult.ratio > 0.99)
        cellXid = [cellXid; j1];
    end
    
end

  CellDif.neuronMeanSpike = data_all;
  
  
  
  
% for j = 1:n
%     for ii = 1:size(neuronBinL,2)-2
%          if neuronBinL(j,ii)  && neuronBinL(j,ii+1) % && neuronBinL(j,ii+2)
%            cellXid = [cellXid; j];
%            break
%            
%        end
%    end
% end
    CellDif.cellXid = cellXid;

end

