% speech enhancement for CHiME4 6-ch data
% corresponding to "The THU-SPMI CHiME-4 system : Lightweight design with
% advanced multi-channel processing, feature enhancement, and language modeling"

% based on CHiME3 official code
% Jon Barker, Ricard Marxer, Emmanuel Vincent, and Shinji Watanabe, The
% third 'CHiME' Speech Separation and Recognition Challenge: Dataset,
% task and baselines, submitted to IEEE 2015 Automatic Speech Recognition
% and Understanding Workshop (ASRU), 2015.

% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

addpath('../enhan');
addpath('../utils');

clear;
micPosition=[0,-0.0385,0;
             0,+0.0385,0;
             -0.033341,-0.01925,0;
             0.033341,0.01925,0;
             -0.033341,+0.01925,0;
             0.033341,-0.01925,0;
             ];
% create template to find the DOA; 
numberLevels = 3;
mySphere=sphereInit(numberLevels);         
Nsource = 2;
Nchan = 6;
Lwindow = 256;
Nfft = Lwindow;
overlap = 0.5;
powThresh = -10;
cmin = 6400; % minimum context duration (400 ms)
cmax = 12800; % maximum context duration (800 ms)
EMITERNUM = 20;
GAMMA = 20;
sampRate=16e3;
num_Fbin = round(Lwindow/2)+1;


% create steering vector template;
waveSpeed=343;
steerVec=zeros(Nchan,num_Fbin,mySphere.SPHERE_NUMBERPOINTS);

%get the steering vector;
for j=1:mySphere.SPHERE_NUMBERPOINTS
    delay = sqrt(sum((mySphere.spherePoints(j,:)-micPosition).*...
    (mySphere.spherePoints(j,:)-micPosition),2))/waveSpeed;%broadcast
    delay=delay-min(delay);
    for index=1:num_Fbin
        f=(index-1)*(sampRate/2)/(Lwindow/2);
        steerVec(:,index,j)=exp(-2*pi*1i*f*delay);
    end    
end

dirn='D:\out16k\midata-time-freq-mask-beamforming\';
filename = 'person2_scene1_music1_level1_angle3_noise9_SNR-5dB_01.wav';
 xsize = size(audioread([dirn filename]));
 x = audioread([dirn filename]);
                  
[ftbin,Nframe,Nbin,Lspeech] =  STFT(x, Lwindow, overlap, Nfft);   
% for GSC fixed beamformer
targetY = squeeze(mean(ftbin,1));

% correlation of X
XX = bsxfun(@times, permute(ftbin,[1,4,2,3]), conj(permute(ftbin,[4,1,2,3])));

Xcor = mean(XX, 4);
% Ncor = mean(XX(:,:,:,[1:10,Nframe-10:Nframe]),4);
% Load context (up to 5 s immediately preceding the utterance)
% noise = read_context(nchan, c_ind, Path, mode, mat, real_mat, uttInd, cmax, cmin);            
softmask = cGaussMask(ftbin,Nsource,XX,EMITERNUM);
Ncor = bsxfun(@rdivide, mean(bsxfun(@times, XX, ...
    permute(softmask(:,:,2),[3,4,1,2])),4),permute(mean(softmask(:,:,2),2), [2,3,1]));                       
Gcor = Xcor - Ncor;


[Y,Df] = MVDR_EV(ftbin, Gcor, Ncor);

%
%match_values of SpherePoints; 
Pmatch=zeros(mySphere.SPHERE_NUMBERPOINTS,4);
Pmatch(:,1:3) = mySphere.spherePoints;

% template match
for j=1:mySphere.SPHERE_NUMBERPOINTS
        Pmatch(j,4) = sum(sum((Df.*conj(permute(steerVec(:,:,j),[1,2])))).*...
        conj(sum(Df.*conj(permute(steerVec(:,:,j),[1,2])))));
end

% plot the sphere;
figure;scatter3(Pmatch(:,1),Pmatch(:,2),Pmatch(:,3),256,Pmatch(:,4));

%

%delaySum use the Df;
yDelaySumDfFt = delaySum(ftbin,Df);
yDelaySumDf = ISTFT(yDelaySumDfFt,Lwindow,overlap);
yDelaySumDf = yDelaySumDf/max(abs(yDelaySumDf));
%delay Sum use the DOA got by match template
index = find(Pmatch(:,4) == max(Pmatch(:,4)));
bestSteerVec = steerVec(:,:,index);
yDelaySumDoaFt = delaySum(ftbin,bestSteerVec/6);
yDelaySumDoa = real(ISTFT(yDelaySumDoaFt,Lwindow,overlap));
yDelaySumDoa = yDelaySumDoa/max(abs(yDelaySumDoa));
% beamforming by MVDR
yBssMvdrEgFt = MVDR_EV(ftbin, Gcor, Ncor);
yBssMvdrEg = ISTFT(yBssMvdrEgFt,Lwindow,overlap);
yBssMvdrEg = yBssMvdrEg / max(abs(yBssMvdrEg));
% audiowrite([dirn filename  'yBssMvdrEg' '.wav'],yBssMvdrEg, 16000);
figure;subplot(4,1,1);plot(x(:,1)/max(abs(x(:,1))));title('原始信号-离线结果');
subplot(4,1,2);plot(yDelaySumDf);title('DelaySum+steering Vector');
subplot(4,1,3);plot(yDelaySumDoa);title('DelaySum+steering Vector match DOA');
subplot(4,1,4);plot(yBssMvdrEg);title('MVDR+steering Vector');

audiowrite([dirn filename  'x' '.wav'],x(:,1)/max(abs(x(:,1))), 16000);
audiowrite([dirn filename  'yDelaySumDf' '.wav'],x(:,1)/max(abs(x(:,1))), 16000);
audiowrite([dirn filename  'yDelaySumDf' '.wav'],yDelaySumDf, 16000);
audiowrite([dirn filename  'yDelaySumDoa' '.wav'],yDelaySumDoa, 16000);
audiowrite([dirn filename  'yBssMvdrEg' '.wav'],yBssMvdrEg, 16000);
 