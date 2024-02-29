%% Analysis of spatial information and exploration-dependent place cells
% ZhouLab, iHuman, ShanghaiTech, 07/30/2021
clear all
close all
clc
 
addpath(genpath('helpers'))
addpath(genpath('Example_data'))

load('Data_MS#6.mat')

%% Pre-processing behavior videos to generate anmail moving trajectories
% addpath(genpath('ms_Scripts')) % Code from https://github.com/etterguillaume/MiniscopeAnalysis
 
analyse_behaviory = false;
if analyse_behaviory
    behav = msGenerateVideoObj(pwd,'behavCam');
    behav = msSelectPropsForTracking(behav);
    trackLength = 35; %cm
    behav = msExtractBehavoir(behav, trackLength);
    save('behav_new.mat','behav','-v7.3');

end


%% Generate coordinates of animals trajectory


NeuronP = spike_prob*30;
NeuronP(isnan(NeuronP))=0;
[n, m] = size(NeuronP);
lapNum = size(Lap,1);

CaCam = 0; BehavCam = 1;
tsdata = importdata('timestamp.dat');
traj1 = generate_trajecotry1(behav,tsdata,Lap,CaCam,BehavCam);



%% Find coordinates for reward

  
thetaPosition = traj1.thetaPosition;
   figure
   Post_lapEnd = zeros(size(Lap));
   for ii = 1:numel(Post_lapEnd)
       Post_lapEnd(ii) = thetaPosition(Lap(ii));
       polarscatter(Post_lapEnd(ii) ,1),hold on
   end

thetaR = mean(Post_lapEnd,'all');
polarscatter(thetaR ,1,'+','r'),hold on
hold off
title('Reward Position')

%% Define Bin activity

nbin = 24;

binActivity = difine_bin_activity(NeuronP,traj1,Lap,nbin,thetaR);


%% Pick starting/end points for each exploration period

pick_points_1  % Click on the location of the objects in sequence


%% Identify laps with distinct exploration activity

Obj4analysis = 1;   % Location of the object according to the picking sequence

ObjBin = Crd.pointBin(Obj4analysis,1):Crd.pointBin(Obj4analysis,2);
obj_point = wrapToPi((Crd.pointBin(Obj4analysis,:))*pi*2/nbin);


%% Load the  manually annotated behavior labels

% load('explap.mat')

%% Check the assigned laps with figures

ShowFigure = true;   % True if needing to check trajectories
ShowVideo = false;   % True if needing to check original videos. Videos must be coppied into the current folder.
ShowCell = false;   % True if showing single cell activity

lap4analysis = 3   ;  % Input a lap number
    
if ShowFigure
    classify_laps_showfig1;
end



%% Identify oePCs - fist step: mean difference with bootstrapping methods


CellDif = specify_cell_bs1(NeuronP,explap,Lap);
cellX1a = CellDif.cellXid;
maxSR1 = max(CellDif.neuronMeanSpike,[],2);  % second step: Max activity >= 2
cellX1 = intersect(find(maxSR1>=2),cellX1a);


%% Show trajectory of oePCs - spike map of individual cells

cellID =  cellX1(4) ;


object_spike_map3


%%  Show hotmap of oePCs

cellid =  cellX1(4) ;

lapseq = [explap.explap1; explap.explap2];

figure
clim = [0 0.5];

LapBin = activity_map(binActivity,cellid,lapseq);
imagesc(LapBin)
yline(length(explap.explap1)+0.5,'w')
   
sgtitle(['Cell # ',...
    num2str(cellid),'; exploration laps, 1-', num2str(numel(explap.explap1))])

bsResult1 = bootstrap_analysis(CellDif.neuronMeanSpike(cellid,:),length(explap.explap1),10000);
figure
histogram(bsResult1.mean_b,'Normalization','probability','FaceColor','w'),hold on
xline(bsResult1.mean_ori,'Color','r','LineWidth',2)
xlabel('Activity difference')
ylabel('Fraction of shuffles')
title(['P = ' num2str(1-bsResult1.ratio)])
%% Show SFP of oePCs

% load('NeuronA'); 
% NeuronA = reshape(neuron.A,Ysiz(1),Ysiz(2),[]);  % Data output from CNMF-E

neuronAbg = (max(NeuronA,[],3));

figure

neuronAeg = max(NeuronA(:,:,cellid),[],3);

neuronShow = cat(3,neuronAbg+neuronAeg*3,neuronAbg, neuronAbg);

imshow(neuronShow)
   

sgtitle(['Cell# ' num2str(cellid)])




%% Identification of place cells from all laps

disp('Computing spatial information...')

nbin2 = 24;

binActivity_pc = difine_bin_activity(NeuronP,traj1,Lap,nbin2,thetaR);

SI_original = spatial_information(binActivity_pc); % Calculate original spatial information

% Calculate shuffled spatial information
repeatNum = 100;
SI_shuffle = zeros(n,repeatNum);
for j = 1:repeatNum
    NeuronSL_shf = randshift(binActivity_pc.neuronSL_inlap);
    SI_temp = spatial_information(binActivity_pc, NeuronSL_shf);
    SI_shuffle(:,j) = SI_temp.SIneuron;
end

% To determine if a neuron is a place cell
placecell = significance_determine(SI_original.SIneuron,SI_shuffle,2.3);

PlaceCell_lg = placecell.Determinant & SI_original.neuronOverallRate > mean(SI_original.neuronOverallRate)* 0.1; 
PlaceCell_id = find(PlaceCell_lg==1);
PlaceCell_ratio = sum(PlaceCell_lg)/numel(PlaceCell_lg)

activityMapAll = SI_original.neuronBinRate; 


%% Compute place cells for explored or unexplored laps only

disp('Computing spatial information...')

repeatNum = 100;
binActivity_exp = difine_bin_activity(NeuronP,traj1,Lap(explap.explap1,:),nbin2,thetaR);
SI_exp = spatial_information(binActivity_exp);
SI_shuffle = zeros(n,repeatNum);
for j = 1:repeatNum
    NeuronSL_shf = randshift(binActivity_exp.neuronSL_inlap);
    SI_temp = spatial_information(binActivity_exp, NeuronSL_shf);
    SI_shuffle(:,j) = SI_temp.SIneuron;
end
placecell_exp = significance_determine(SI_exp.SIneuron,SI_shuffle,2.33);
SSI_exp = placecell_exp.ZValue;
activityMap_exp = SI_exp.neuronBinRate; 



binActivity_unexp = difine_bin_activity(NeuronP,traj1,Lap(explap.explap2,:),nbin2,thetaR);
SI_unexp = spatial_information(binActivity_unexp);
SI_shuffle = zeros(n,repeatNum);
for j = 1:repeatNum
    NeuronSL_shf = randshift(binActivity_unexp.neuronSL_inlap);
    SI_temp = spatial_information(binActivity_unexp, NeuronSL_shf);
    SI_shuffle(:,j) = SI_temp.SIneuron;
end
placecell_unexp = significance_determine(SI_unexp.SIneuron,SI_shuffle,2.33);
SSI_unexp = placecell_unexp.ZValue;
activityMap_unexp = SI_unexp.neuronBinRate; 

%% Show hotmap of all place cells

figure
for j1 = 1:length(PlaceCell_id)
   subplot(ceil(sqrt(length(PlaceCell_id))),ceil(sqrt(length(PlaceCell_id))),j1)
   acttemp = permute(binActivity_pc.neuronBinLapS(PlaceCell_id(j1),:,:),[3,2,1]);
   imagesc(acttemp)
   title(PlaceCell_id(j1))
end

%% To determine whether cellX are place cells within exploration laps
% Identify oePCss - third step

cellX = cellX1(ismember(cellX1,placecell_exp.ID));

disp(['The identified oePCs are: ', num2str(cellX')])



%% Activity maps with peak bin sequence

activitymap1 = activityMapAll(PlaceCell_id,:);
% activitymap1 = activityMapAll;

activitymap2 = activity_map2(activitymap1);
figure

imagesc(activitymap2.map)

sgtitle('Activity maps of all place cells')
activitymap2.cellID = PlaceCell_id(activitymap2.idx);

%% Center of mass
COM=zeros(n,1);
for j=1:n
   
   BinLap=binActivity.neuronBinLapS(j,:,:);
   LapBin=transpose(reshape(BinLap,24,[]));
   LapBin_aver=mean(LapBin(:,:));
   [C_theta, C_bin ] = center_of_mass(LapBin_aver,nbin);
   COM(j,1)=C_theta;
   COM(j,2)=C_bin;
  
end
COM1 = COM(:,1)+pi;  % range in 0 - 2*pi


%%  Difference index 


bin_PC = Crd.objBin(Obj4analysis)-1:Crd.objBin(Obj4analysis)+1;

explored = activityMap_exp(:,bin_PC);
unexplored = activityMap_unexp(:,bin_PC);

diff_index = zeros(n,1);
    for j = 1:n
        diff_index(j,1)=(max(explored(j,:))-max(unexplored(j,:)))./(max(explored(j,:))+max(unexplored(j,:)));
    end


