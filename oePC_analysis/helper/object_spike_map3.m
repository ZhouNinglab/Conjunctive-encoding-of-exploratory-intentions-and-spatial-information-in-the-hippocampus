%% Show spiking activity of a neuron during exploration of  specific object in all laps
% opt = 1; % show spiking during exploration only
% opt ==2, %show spiking during the entire lap
opt = 3; % show timeing of spikesc1factor = 0.8;
c1factor = 0.8;

x_r = traj1.x1;
y_r = traj1.y1;


t2 = ["Explored laps" "Unexplored laps"];


for i1 =  1:2 
    
if i1 ==1
    Lap1 = [explap.explap1];
elseif i1 ==2
    Lap1 = [explap.explap2]; 
elseif i1 ==3
    Lap1 = lap_misplace; 
end
    
figure ('Position',[100 900-i1*250 500 150])    
    
Lap1 = sortrows(Lap1);
laplabel = {'Exploration','Non-Exploration'};

for j2= 1: size(Lap1,1)

    subplot(1,length(Lap1), j2)
    lapview.lapfms_temp = traj1.Lap(Lap1(j2),1):traj1.Lap(Lap1(j2),2);

    lapview.theta_temp = binActivity.laps_theta_post{Lap1(j2)};
    fms1 = lapview.lapfms_temp(find(lapview.theta_temp>(obj_point(1)-0.1) & lapview.theta_temp< (obj_point(2))));
    plot(x_r(lapview.lapfms_temp),y_r(lapview.lapfms_temp),':k' ,'LineWidth',1),hold on

        plot(x_r(fms1),y_r(fms1),'k' ,'LineWidth',1),axis ij, hold on

     xlim([-50 50])
     ylim([-50 50])
     xticklabels([]);
     yticklabels([]);
c1 = max(binActivity.neuronSL_inlap(cellID,:)) * c1factor;      
    if opt == 1

            C1 =NeuronP(cellID,fms1);

             scatter((x_r(fms1)),(y_r(fms1)),5,C1,'filled'), hold on
             title((Lap1(j2)))

    end
   
    caxis([0 c1])
% colormap jet
    
    if opt ==2
        C1 =NeuronP(cellID,lapview.lapfms_temp);
        scatter(x_r(lapview.lapfms_temp),y_r(lapview.lapfms_temp),15,'k'), hold on
        scatter(x_r(lapview.lapfms_temp),y_r(lapview.lapfms_temp),10,C1,'filled')
        title((Lap1(j2)))
        caxis([0 c1])

    end
    
    if opt == 3
        NeuronSL = logical(NeuronS >= 0.1);
        for i2 = 1:length(fms1)
            if NeuronSL(cellID,fms1(i2))
               scatter(x_r(fms1(i2)),y_r(fms1(i2)),10,'filled','MarkerFaceColor','r')
            end
        end

    end


end
sgtitle(laplabel{i1})

end