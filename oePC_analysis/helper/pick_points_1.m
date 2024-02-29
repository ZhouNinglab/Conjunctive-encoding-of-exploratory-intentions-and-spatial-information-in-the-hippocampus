%% Pick the starting and ending points of object exploration

figure


plot(traj1.x1,traj1.y1),axis ij, hold on
[Crd.xR, Crd.yR] = pol2cart(thetaR,median(traj1.rhoPosition,'omitnan'));
scatter(Crd.xR, Crd.yR, 100,'+','r')

disp('Click to pick the locations of objects, double click to finish.')
[Crd.x2s,Crd.y2s] = getpts;
scatter(Crd.x2s,Crd.y2s,50,'r','filled')



%% Difine bins for each point
Crd.xys = [Crd.x2s,Crd.y2s];
% Crd.xys = [Crd.x2s,Crd.y2s; Crd.x2e,Crd.y2e; Crd.x3s,Crd.y3s; Crd.x3e,Crd.y3e];
[Crd.objtheta, Crd.objrho] = cart2pol(Crd.xys(:,1),(Crd.xys(:,2)));
Crd.theta_post = zeros(size(Crd.objtheta));
Crd.theta_pre = zeros(size(Crd.objtheta));



Crd.theta_pre = thetaR - Crd.objtheta;
Crd.theta_post = wrapTo2Pi(Crd.theta_pre);
Crd.edge1 = 0:2*pi/nbin:2*pi;
Crd.objBin = discretize(Crd.theta_post,Crd.edge1);

Crd.pointBin = zeros(length(Crd.theta_post),2);
Crd.pointBin(:,1) = Crd.objBin - 3;
Crd.pointBin(:,2) = Crd.objBin + 3;

