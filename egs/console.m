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

clear all;

Nsource = 2;
Nchan = 6;
Lwindow = 256;
Nfft = Lwindow;
overlap = 0.75;
powThresh = -10;
cmin = 6400; % minimum context duration (400 ms)
cmax = 12800; % maximum context duration (800 ms)
EMITERNUM = 20;
GAMMA = 20;

dirn = 'C:\Users\mi\Downloads\myairbridge-FDXXE174Y3RM9Q25\CHiME3_isolated_et05_real\CHiME3\data\audio\16kHz\isolated\et05_bus_real\';
firname= 'F05_440C020E_BUS';
% Load WAV files
xsize = size(audioread([dirn firname '.CH1.wav']));
x = zeros(xsize(1),Nchan);
xc = zeros(xsize(1),Nchan);
for clp = 1:Nchan,
    [x(:,clp),fs] = audioread([dirn firname '.CH' int2str(clp) '.wav']);
    xc(:,clp) = x(:,clp) / norm(x(:,clp)) * norm(x(:,1));
end

% Check microphone failure
CLwindow = 2048;
Coverlap = 0.5;
CNfft = CLwindow;
[~,~,~,~,CspeechFrame] =  STFT(xc, CLwindow, Coverlap, CNfft);
CspeechFrame = permute(CspeechFrame, [3,1,2]);
NcorFFT = 2 * (CNfft-1);
c1 = fft(repmat(permute(CspeechFrame, [1,4,2,3]),[1 Nchan 1 1]), NcorFFT, 3);
c2 = fft(repmat(permute(CspeechFrame, [4,1,2,3]),[Nchan 1 1 1]), NcorFFT, 3);
corr = real(ifft(c1 .* conj(c2),NcorFFT,3));
corr = squeeze(max(corr, [], 3));

corrjudge = sum(corr, 2) - sum(bsxfun(@times,corr,eye(Nchan)),2);
corrjudge = squeeze(corrjudge);
corrjudge_rate = bsxfun(@rdivide, corrjudge, median(corrjudge,1));
fail = corrjudge_rate < 0.6;
fail = (sum(fail,2) > 1);
fail = fail';
if all(fail)
    fail = (sum(corrjudge_rate,2) == max(sum(corrjudge_rate,2)));
    fail = ~fail;
    fail = fail';
end
% choose reference mic
[~,refMic] = max(sum(corrjudge(~fail,:),2));          

Nchan = sum(~fail);
x = x(:,~fail);
            
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
Ncor = bsxfun(@rdivide, mean(bsxfun(@times, XX, permute(softmask(:,:,2),[3,4,1,2])),4),permute(mean(softmask(:,:,2),2), [2,3,1]));                       
Gcor = Xcor - Ncor;
