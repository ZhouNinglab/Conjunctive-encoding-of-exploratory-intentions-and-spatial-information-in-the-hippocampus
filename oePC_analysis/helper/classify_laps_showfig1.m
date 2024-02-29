%% Characterization of behavioural bouts


x_r = traj1.x1;
y_r = traj1.y1;
j2 = lap4analysis;
lapview.lapfms_temp = traj1.Lap(j2,1):traj1.Lap(j2,2);
lapview.theta_temp = binActivity.laps_theta_post{j2};
fms1 = lapview.lapfms_temp(find(lapview.theta_temp>obj_point(1)-0.2 & lapview.theta_temp<obj_point(2)));

    
lapview.rho_temp = traj1.rhoPosition(fms1);
lapview.speed_temp0 = transpose(smooth(traj1.speed));
lapview.speed_temp = lapview.speed_temp0(fms1);
% lapview.speed_temp = (behav.speed1(fms1))';
lapview.drho = diff(lapview.rho_temp);
lapview.dspeed = diff(lapview.speed_temp);

lapview.fms1 = fms1;
lapview.points = lapview.drho>0.15 & lapview.dspeed <-1;
lapview.fms2 = fms1(lapview.points);
sum(lapview.points);
%% Trajectory of behavioral bout

figure('position', [200 200 700 700])
sgtitle(['Lap ',num2str(j2), '; Position ', num2str(Obj4analysis)])


subplot(331)
plot(x_r(lapview.lapfms_temp),y_r(lapview.lapfms_temp),':k' ),hold on
Check_trajectory
hold on

scatter(x_r(lapview.fms2),y_r(lapview.fms2),10,'filled','MarkerFaceColor','g')
title('Behavioral bout')
%% Show Speed & Rho over frames
subplot(332)

 plot(lapview.speed_temp),hold on
plot(lapview.rho_temp),hold off
legend('Speed','Rho')
title('Speed & Rho / frames')


%% Show Speed with trajectory
 subplot(333)
 plot(x_r(lapview.lapfms_temp),y_r(lapview.lapfms_temp),':k' ),hold on
plot(x_r(lapview.fms1),y_r(lapview.fms1),'k')
scatter(x_r(lapview.fms1),y_r(lapview.fms1),10,1:length(lapview.fms1),'filled')
axis ij
xlim([-50 50])
ylim([-50 50])
colormap copper
freezeColors;
title('Speed with trajectory')

%% Speed/Time relationship in each bout
subplot(334)
% fms1 = explap.fms{lapseq(j2)};
        rho_temp = traj1.rhoPosition(fms1);
        speed_temp = traj1.speed(fms1)/(behav.trackLength/35);
        theta_temp = traj1.thetaPosition(fms1);
        edgetemp = (floor(min(10*theta_temp)):1:ceil(max(10*theta_temp)))/10;
        tempbin = discretize(theta_temp,edgetemp);
        lap_speed = zeros(length(edgetemp)-1,1);
        lap_time = zeros(length(edgetemp)-1,1);
        for j3 = 1:length(edgetemp)-1
           tempid = find(tempbin==j3);
           lap_speed(j3) = mean(speed_temp(tempid));
           lap_time(j3) = length(tempid)/30;
        end

        plot(lap_speed,lap_time,'k'), hold on
        scatter(lap_speed,lap_time,20,(1:length(edgetemp)-1)','filled'),colormap copper
        freezeColors;
        xlim([0 120])
        ylim([0 1.5])
        title('Speed/Time')
ylim([0 inf])
xlim([0 inf])
%% Check with videos
if ShowVideo
     lapview.vn1 = ceil(fms1(1) / 1000);
     lapview.vn2 = ceil(fms1(end) / 1000);
     lapview.ftemp = mod(fms1,1000); 
     if lapview.ftemp(end) == 0
         lapview.ftemp(end) =  1000;
     end
     
         if lapview.vn1 == lapview.vn2
             lapview.VideoFile1 = ['behavCam' num2str(lapview.vn1) '.avi'];
             lapview.v1 = VideoReader(lapview.VideoFile1);
             lapview.video1 = read(lapview.v1,[lapview.ftemp(1) lapview.ftemp(end)]);


         else

             lapview.VideoFile1a = ['behavCam' num2str(lapview.vn1) '.avi'];
             lapview.v1a = VideoReader(lapview.VideoFile1a);
             lapview.video1a = read(lapview.v1a,[lapview.ftemp(1) 1000]);
             lapview.VideoFile1b = ['behavCam' num2str(lapview.vn2) '.avi'];
             lapview.v1b = VideoReader(lapview.VideoFile1b);
             lapview.video1b = read(lapview.v1b,[1 lapview.ftemp(end)]);

             lapview.video1 = cat(4, lapview.video1a, lapview.video1b);

         end
    
          lapview.v2 = VideoWriter('Trajectory4test.avi');
          lapview.v2.FrameRate = 20;
          open(lapview.v2) % videowriter file needs to be opened before write
          writeVideo(lapview.v2,lapview.video1) %write the data from the trancated video into the videowriter file
          close(lapview.v2)

         implay('Trajectory4test.avi')
end

%% Calculate off-track points


YL = [25 50];

fmsexplap = cell(1,2);
rhoexplap = cell(1,2);
thetaexplap = cell(1,2);
for i1 = 1:2
        if i1 == 1
         laptemp = lap4analysis;
        elseif i1 == 2
             laptemp =  explap.explap2;
        end
        
%          rhotemp = [];
    fmstemp = [];           
%     thetatemp = [];
    for j2=1:length(laptemp)
         j1 = laptemp(j2);
    %      rhotemp = [rhotemp traj1.rhoPosition(Lap(j1,1):Lap(j1,2))];
         fmstemp = [fmstemp explap.fms{j1}];
    %      fmstemp = [fmstemp Lap(j1,1):Lap(j1,2)];
    %    plot(traj1.x1(Lap(j1,1):Lap(j1,2)),traj1.y1(Lap(j1,1):Lap(j1,2)),'b'),hold on 
    end
fmsexplap{i1} = fmstemp;
rhoexplap{i1} = traj1.rhoPosition(fmsexplap{i1});
thetaexplap{i1} = traj1.thetaPosition(fmsexplap{i1});
end




edge1 = (floor(min(10*horzcat(thetaexplap{:,:}))):1:ceil(max(10*horzcat(thetaexplap{:,:}))))/10;
expbin = cell(1,2);
expbin{1,1} = discretize(thetaexplap{1},edge1);
expbin{1,2} = discretize(thetaexplap{2},edge1);
upperlim = zeros(1,length(edge1)-1);

tempid = [];
for j2 = 1 :length(edge1)-1
    temp1 = find(expbin{1,2}==j2);
    temp2 = rhoexplap{1,2};
    upperlim(j2) = mean(temp2(temp1)) + 2 * std(temp2(temp1));
    temp4 = rhoexplap{1,1};
    temp3 = find(expbin{1,1}==j2);
    tempid = [tempid temp3(temp4(temp3)<=upperlim(j2))];
    
end
tempid = transpose(sortrows(tempid'));
temp = fmsexplap{1,1};
fmsexplap_n = temp(tempid);
temp = rhoexplap{1,1};
rhoexplap_n = temp(tempid);
temp = thetaexplap{1,1};
thetaexplap_n = temp(tempid);



     
%% Show Trajectories of expl. bout (red) with non-expl. bouts (blue)
 subplot(335)
 for j1 = 1:length(explap.explap2)
     laptemp = explap.explap2(j1);
      plot(traj1.x1(Lap(laptemp,1):Lap(laptemp,2)),traj1.y1(Lap(laptemp,1):Lap(laptemp,2)),'-k'),hold on
 end
 
    plot(traj1.x1(fmsexplap{1,2}),traj1.y1(fmsexplap{1,2}),'b'),hold on
    plot(traj1.x1(Lap((lap4analysis),1):Lap((lap4analysis),2)),traj1.y1(Lap((lap4analysis),1):Lap((lap4analysis),2)),'-.r')
    plot(traj1.x1(fmsexplap{1,1}),traj1.y1(fmsexplap{1,1}),'r'),hold on
    axis ij
    
    if exist('Fms_nearobj')==1
    hold on, scatter(traj1.x1(Fms_nearobj(lap4analysis)),traj1.y1(Fms_nearobj(lap4analysis))) % show near object position
    end
    title('Trajectories of expl. bout')

%% Show 'Off-track points'
subplot(336)
      plot(thetaexplap{1,1},rhoexplap{1,1},'r'),hold on
      scatter(thetaexplap{1,1},rhoexplap{1,1},20,'r'),hold on
     scatter(thetaexplap_n,rhoexplap_n,20,'r','filled')
     plot(thetaexplap{1,2},rhoexplap{1,2},'b')
%      scatter(thetaexplap{1,2},rhoexplap{1,2},20,'b','filled','^')
     scatter(thetaexplap{1,2},rhoexplap{1,2},20,'b','filled','^')
     plot(edge1(2:end)-0.05,upperlim,'--k','LineWidth',2)
     set(gca,'XDir','reverse')
     ylim(YL)
     title('Off-track points')
     
     
%% Show Speed in theta-rho plot
subplot(337)
plot(thetaexplap{1,1},rhoexplap{1,1},'-k'),hold on
     scatter(thetaexplap{1,2},rhoexplap{1,2},20,traj1.speed(fmsexplap{1,2})/(behav.trackLength/35),'filled','^')
%      scatter(thetaexplap_n,rhoexplap_n,20,traj1.speed(fmsexplap_n),'filled'), hold on
     scatter(thetaexplap{1,1},rhoexplap{1,1},20,traj1.speed(fmsexplap{1,1})/(behav.trackLength/35),'filled')
     plot(edge1(2:end)-0.05,upperlim,'--k')
     set(gca,'XDir','reverse')  
     colorbar
     colormap pink
     freezeColors; freezeColors(colorbar)
     title('Speed')
      ylim(YL)

%% Show cell activity with trajectory
if ShowCell 
   id8 =  cellid;

         
      subplot(338)
      plot(thetaexplap{1,1},rhoexplap{1,1},'-k'),hold on
     scatter(thetaexplap{1,2},rhoexplap{1,2},20,NeuronP((id8),fmsexplap{1,2}),'filled','^')
%      scatter(thetaexplap_n,rhoexplap_n,20,NeuronP(cellX1(id8),fmsexplap_n),'filled'), hold on
     scatter(thetaexplap{1,1},rhoexplap{1,1},20,NeuronP((id8),fmsexplap{1,1}),'filled'), hold on
     
     plot(edge1(2:end)-0.05,upperlim,'-k')
     set(gca,'XDir','reverse')  
     ylim(YL)
     colormap default

     CL = [0 quantile([NeuronP((id8),fmsexplap{1,1}) NeuronP((id8),fmsexplap{1,2})],0.90)];
     caxis(CL)
%     colorbar
     freezeColors; freezeColors(colorbar)
     title(id8)
     
    
    subplot(339)
%       scatter(fmsexplap{1,1},rhoexplap{1,1},20,NeuronP((id8),fmsexplap{1,1}),'filled','LineWidth',1,'MarkerEdgeColor','k'), hold on
plot(fmsexplap{1,1},traj1.speed(fmsexplap{1,1})), hold on
      scatter(fmsexplap{1,1},traj1.speed(fmsexplap{1,1}),20,NeuronP((id8),fmsexplap{1,1}),'LineWidth',1) %,'filled','LineWidth',1,'MarkerEdgeColor','k'), hold on
colormap default
% plot(fmsexplap_n,traj1.speed(fmsexplap_n))
scatter(fmsexplap_n,traj1.speed(fmsexplap_n),20,NeuronP((id8),fmsexplap_n),'filled') %,'LineWidth',1,'MarkerEdgeColor','w')
 title('Inferred activity')
caxis(CL)
  colorbar
  freezeColors; freezeColors(colorbar)
  



end
