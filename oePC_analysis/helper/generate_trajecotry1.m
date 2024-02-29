function [trajectory] = generate_trajecotry(behav,tsdata,Lap,Cam0,Cam1)
%GENERATE_TRAJECOTRY 
%   produce time-dependent coordinates of animals movement




data1 = tsdata.data;
ts0 = data1(data1(:,1)==Cam0,:);
ts1 = data1(data1(:,1)==Cam1,:);
ts0(1,3) = 1;
ts1(1,3) = 1;
m = size(ts0,1);
% m1 = size(ts1,1);
t = ts1(:,3);
timeq = ts0(:,3);

while length(unique(t)) ~= length(t)
    [~,ic,~] = unique(t);
    t(~ismember((1:length(t)),ic),:) = t(~ismember((1:length(t)),ic),:)+1;
end

if length(unique(timeq)) ~= length(timeq)
    [~,ic,~] = unique(timeq);
    timeq(~ismember((1:length(timeq)),ic),:) = timeq(~ismember((1:length(timeq)),ic),:)+1;
end


v1 = behav.position(:,1);
v2 = behav.position(:,2);
xx = interp1(t,v1,timeq);
yy = interp1(t,v2,timeq);
speedraw = interp1(t,behav.speed,timeq);
Laptemp = reshape(Lap,[],1);
Ttemp1 = interp1(ts1(:,2),t,Laptemp);
Lapnew = interp1(timeq,ts0(:,2),Ttemp1);
Lapnew = reshape(round(Lapnew),[],2);
% 
% xx = interp1(positionraw(:,1),(1:m));
% yy = interp1(positionraw(:,2),(1:m));
% speedraw = interp1(speedraw, (1:m));


xCenter=mean([max(xx) min(xx)]);
yCenter=mean([max(yy) min(yy)]);  % find center coordinates of x and y
x1=xx-xCenter;
y1=yy-yCenter;

[thetaPosition,rhoPosition]=cart2pol(x1,y1);   

trajectory.xraw = xx';
trajectory.yraw = yy';

trajectory.x1 = x1';
trajectory.y1 = y1';
trajectory.thetaPosition = thetaPosition';
trajectory.rhoPosition = rhoPosition';
trajectory.speed = speedraw';
trajectory.Lap = Lapnew;
end

