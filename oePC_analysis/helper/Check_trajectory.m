%% Draw moving trajectory
x_r = traj1.x1;
y_r = traj1.y1;
% x_r = traj1.xraw;
% y_r = traj1.yraw;

% figure

 xlim([-50 50])
 ylim([-50 50])
axis ij
h = animatedline;
for ii= fms1; 
    addpoints(h,x_r(ii),y_r(ii));

    drawnow % limitrate
end
