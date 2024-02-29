function [traj2,NeuronSL2,Lap2] = speed_filter(traj1,NeuronSL,Lap,spdthr)
%SPEED_FILTER 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin < 4
    spdthr = 10;
end
speed1 = traj1.speed;
fms = (speed1 >= spdthr);
x1 = traj1.x1;
y1 = traj1.y1;
thetaPosition = traj1.thetaPosition;
rhoPosition = traj1.rhoPosition;
traj2.x1 = x1(fms);
traj2.y1 = y1(fms);
traj2.thetaPosition = thetaPosition(fms);
traj2.rhoPosition = rhoPosition(fms);
traj2.speed = speed1(fms);

lapIdx = zeros(size(speed1));
lapNum = size(Lap,1);
for j = 1:lapNum
   lapIdx(Lap(j,1):Lap(j,2))=j;
end
 
lapIdx = lapIdx(fms);
Lap2 = zeros(max(lapIdx),2);
for ii = 1:max(lapIdx)
    Lap2(ii,1) = min(find(lapIdx == ii));
    Lap2(ii,2) = max(find(lapIdx == ii));
end

NeuronSL2 = NeuronSL(:,fms);

end

