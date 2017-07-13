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

dirn='./data/';
filename = 'person1_scene1_music1_level1_angle2_noise7_SNR15dB.wav';
 xsize = size(audioread([dirn filename]));
 data = audioread([dirn filename]);
 
 % segment data to get DOA by mask
 len_Batch = 1024*4;
 batch_OverLap = 0.5;
 number_Lwin = round(xsize(1)/(Lwindow*(1-overlap)));
 angleSet = zeros(number_Lwin-30,1);
 for Index = 1:number_Lwin-30
         try
            x = data(1+(Index-1)*Lwindow*(1-overlap)...
                 :len_Batch+(Index-1)*Lwindow*(1-overlap),:);
        [ftbin,Nframe,Nbin,Lspeech] =  STFT(x, Lwindow, overlap, Nfft);   

        % correlation of X
        XX = bsxfun(@times, permute(ftbin,[1,4,2,3]), conj(permute(ftbin,[4,1,2,3])));

        Xcor = mean(XX, 4);
        softmask = cGaussMask(ftbin,Nsource,XX,EMITERNUM);
        Ncor = bsxfun(@rdivide, mean(bsxfun(@times, XX, ...
            permute(softmask(:,:,2),[3,4,1,2])),4),permute(mean(softmask(:,:,2),2), [2,3,1]));                       
        Gcor = Xcor - Ncor;

        %match_values of SpherePoints; 
        Pmatch=zeros(mySphere.SPHERE_NUMBERPOINTS,4);
        Pmatch(:,1:3) = mySphere.spherePoints;

        % get the steering vector which is estimated by time-freq mask,
        % meanwhile the time-freq mask is estimated by CGMM;
        [Y,Df] = MVDR_EV(ftbin, Gcor, Ncor);

        % template match
        for j=1:mySphere.SPHERE_NUMBERPOINTS
                Pmatch(j,4) = sum(sum((Df.*conj(permute(steerVec(:,:,j),[1,2])))).*...
                conj(sum(Df.*conj(permute(steerVec(:,:,j),[1,2])))));
        end

        % plot the sphere;
        % scatter3(Pmatch(:,1),Pmatch(:,2),Pmatch(:,3),256,Pmatch(:,4));

        %find the DOA from Spherepoints
        index = find(Pmatch(:,4) == max(Pmatch(:,4)));
        bestPoint = mySphere.spherePoints(index,:);
        angle = atan(abs(bestPoint(2)/bestPoint(1)))*180/pi;
        if bestPoint(1) < 0 &&bestPoint(2) < 0
            angle = -angle-90;
        end
        if bestPoint(1) < 0 &&bestPoint(2) > 0
            angle = angle+90;
        end
        if bestPoint(1) > 0 &&bestPoint(2) < 0
            angle = -angle;
        end
        angleSet(Index) = angle;
%     fprintf('angle:%f\n',angle);
         end
     continue;
 end



