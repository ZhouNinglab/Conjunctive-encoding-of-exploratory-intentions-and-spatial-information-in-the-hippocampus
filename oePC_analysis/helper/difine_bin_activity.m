function [binActivity] = difine_bin_activity(NeuronSL,traj,Lap,binNum,thetaR)
%DIFINE_BIN_ACTIVITY 
%   Default bin number: 48
%   Default starting postition(theta) 0.8634
% if nargin <=4
%     thetaR = 0.8634;
%     if nargin <=3
%         binNum = 48;
%     end
% end

if floor(binNum)~=ceil(binNum)
    disp('Bin number must be an integer!')
    return
end

[n,m] = size(NeuronSL);
[lapnum,~]=size(Lap);
nbin = binNum;
binwidth = 2*pi/nbin;
% neuronBinLapS = zeros(n, nbin, lapnum);
% speedBinLap = zeros(nbin, lapnum);
% rhoBinLap =  zeros(nbin, lapnum);
speedBinLap = NaN(nbin, lapnum);
rhoBinLap =  NaN(nbin, lapnum);
binvalueAllc = cell(1,lapnum);
NeuronS_allLap = [];
thetaPosition = traj.thetaPosition;
speedraw = traj.speed;
laps_theta_post = cell(lapnum,1);
for jj = 1:lapnum
    theta_thislap = thetaPosition(Lap(jj,1):Lap(jj,2));
    NeuronS_thislap = NeuronSL(:,Lap(jj,1):Lap(jj,2));
    theta_pre = thetaR - theta_thislap;
    theta_post = theta_pre;
    for ii=1:length(theta_thislap)
        if theta_pre(ii)<0
            theta_post(ii)=theta_pre(ii)+2*pi;
        elseif theta_pre(ii)>2*pi
            theta_post(ii)=theta_pre(ii)-2*pi;
        end
    end
    
    laps_theta_post{jj} = theta_post;
    
    binvalue = zeros(length(theta_thislap),1);
    
    speed_thislap = speedraw(Lap(jj,1):Lap(jj,2));

    rho_thislap = traj.rhoPosition(Lap(jj,1):Lap(jj,2));
           
    for j=1:nbin
       
       binvalue(find(theta_post>=(j-1)*binwidth & theta_post<j*binwidth))=j; 
        if ismember(j,binvalue)

       neuronBinLapS(:,j,jj) = mean(NeuronS_thislap(:,find(binvalue==j)),2) ;
%        neuronBinLapS(:,j,jj) = max(NeuronS_thislap(:,find(binvalue==j)),[],2) ;

        speedBinLap(j,jj) = mean(speed_thislap(find(binvalue==j))');  % min
        rhoBinLap(j,jj) = max(rho_thislap(find(binvalue==j))');
        else continue
        end
    end
%     binvalueAll=[binvalueAll binvalue'];
     binvalueAllc{jj} = binvalue';
    neuronBinLapS(isnan(neuronBinLapS))=0;
        
    NeuronS_allLap=[NeuronS_allLap NeuronS_thislap];
end

speedBinLap = speedBinLap';


binActivity.binvalueall = binvalueAllc;
binActivity.binvalue = horzcat(binvalueAllc{:,:});
binActivity.neuronSL_inlap = NeuronS_allLap;
binActivity.speed= speedBinLap;
binActivity.neuronBinLapS = neuronBinLapS;
binActivity.nbin = binNum;
binActivity.rhomax = rhoBinLap';
binActivity.laps_theta_post = laps_theta_post;
end

