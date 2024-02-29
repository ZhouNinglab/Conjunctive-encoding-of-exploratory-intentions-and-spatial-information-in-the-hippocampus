function [SpatialInfo] = spatial_information(binActivity,neuronSL_inlap)
%SPATIAL_INFORMATION 此处显示有关此函数的摘要
%   Default neuronSL_inlap is extracted from binActivity
    if nargin ==1 && isstruct(binActivity)
        neuronSL_inlap = binActivity.neuronSL_inlap;
    end

[n, nbin, lapNum] = size(binActivity.neuronBinLapS);
 binvalueAll = binActivity.binvalue;
 nbin = binActivity.nbin;
 
 P=zeros(1,nbin);
 
for i=1:nbin
     P(i)=length(binvalueAll(binvalueAll==i))/length(binvalueAll);
end

SI=zeros(n,nbin);
neuronBinRate=zeros(n,nbin);

for j=1:n
    for i=1:nbin
        neuronBinRate(:,i)=mean(neuronSL_inlap(:,find(binvalueAll==i)),2);         % n * nbin matrix showing mean NeuronZ in each bin of each neruon
        neuronOverallRate=mean(neuronSL_inlap,2);   % Overall rate of each neuron
        SItemp = P(i) * neuronBinRate(j,i)/ neuronOverallRate(j) * (log2(neuronBinRate(j,i)/ neuronOverallRate(j)));       %SI for each bin
        SItemp(isnan(SItemp))=0;
        SI(j,i)=SItemp;
    end
end

SIneuron  = sum(SI,2);
SpatialInfo.SIneuron = SIneuron;
SpatialInfo.neuronOverallRate = neuronOverallRate;
SpatialInfo.neuronBinRate = neuronBinRate;



end

