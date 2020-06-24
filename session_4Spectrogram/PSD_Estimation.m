% To compute PSD with moving window Hamming, Hanning.
% Author: Rinu Sebastian
% Date: 11/26/2018
% Status: can depend fully
% Ref for normalisation and multiplication with factor 2:
% Algorithm section of https://www.mathworks.com/examples/signal/mw/signal-ex20286697-periodogram
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
% AC Freq: 60 Hz
close all;
%Directory infos
base_dir=pwd;
out_dir = fullfile(base_dir,'figs_PSD'); %creating directory for output if one does not exist
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end

%ecog data and initializations
ecog_struct=load('bp_mot_t_h.mat');
ecog_data=ecog_struct.data; 
clear ecog_struct;
ecog_data= 0.0298*ecog_data; %Convert to microvolt
[time_pts,channel_num]=size(ecog_data); % obtaining time value and # of channels from data
Fs=1000;   % Sampling Frequency
t=(1:time_pts)/Fs; %time vector for time in seconds
N=2^nextpow2(time_pts); %Time point for FFT
ecog_fft=zeros([N,channel_num]);
Pwr_fft=zeros([N/2+1,channel_num]);
pxx=zeros([N/2+1,channel_num]);

%clipping 
Fnotch = 60; % AC Frequency to be notched.
freq= Fs*(0:(N/2))/N; %freq vector for PSD 
clip_freq=find(freq>1*Fnotch-5 & freq<1*Fnotch+5);
for aa=2:round((N/2)/Fnotch)
    temp = find(freq>(aa*Fnotch-5)& freq<(aa*Fnotch+5)); %indices corresponding to 60 Hz or its Harmonics +/- 5 Hz
    clip_freq=[clip_freq temp]; %frequency vector to be clipped.
end 

%with windows
% Hamming/ Hanning window => element wise windowing-> zero padding & FFT calculation
% -> PSD estimation -> clipping -> plotting
over_lap=50; % 50% percentage overlap
seg_length=1000; %dividing signal into segments of 1000, essentially giving 1s segments
nfft=1024; % 1024 point DFT
f=Fs*(0:(nfft/2))/nfft; %Freq vector for Windowed PSD calculation
% clipping
clip_freq_1=find(f>1*Fnotch-5 & f<1*Fnotch+5);
for aa=2:round((nfft/2)/Fnotch)
    temp = find(f>(aa*60-5)& f<(aa*60+5)); %indices corresponding to 60 Hz or its Harmonics +/- 5 Hz
    clip_freq_1=[clip_freq_1 temp]; %frequency vector to be clipped.
end 
pxx_welch_ham = zeros([nfft/2+1,channel_num]);
pxx_welch_han = zeros([nfft/2+1,channel_num]);
[PSD_ham,PSD_han,off_set]=PSD_calc(ecog_data,channel_num,time_pts,over_lap,Fs); %fn call to user defined fn to calculate PSD 
% crosschecking with pwelch()
for kk=1:channel_num
    pxx_welch_ham(:,kk) = pwelch(ecog_data(:,7),hamming(seg_length),off_set,nfft,Fs);
    pxx_welch_han(:,kk) = pwelch(ecog_data(:,7),hanning(seg_length),off_set,nfft,Fs);
end
params.Fs=Fs;
param.pad=1;
params.tapers=[5 9];
[S,freq] = mtspectrumc(ecog_data, params );
for aa=1:round((N/2)/Fnotch)
    temp = find(freq>(aa*Fnotch-5)& freq<(aa*Fnotch+5)); %indices corresponding to 60 Hz or its Harmonics +/- 5 Hz
    clip_freq=[clip_freq temp]; %frequency vector to be clipped.
end 

PSD_ham(clip_freq_1,:)=NaN;
PSD_han(clip_freq_1,:)=NaN;
pxx_welch_ham(clip_freq_1,:)=NaN;
pxx_welch_han(clip_freq_1,:)=NaN;
S(clip_freq,:)=NaN;


h2=figure(3); % plotting windowed PSD
ax2=define_position(2,1);
plot(ax2(2,1),f,10*log10(PSD_ham(:,7)));
ylim(ax2(2,1),[-80 65]);
ylabel(ax2(2,1),{'PSD-Hamming','Power/Freq(dB/Hz)'});
xticklabels(ax2(2,1),[]);
plot(ax2(1,1),f,10*log10(pxx_welch_ham(:,7)));
ylim(ax2(1,1),[-80 65]);
ylabel(ax2(1,1),{'PSD- Hanning','Power/Freq(dB/Hz)'});
xlabel(ax2(1,1),'Frequency (Hz)')
saveas(h2,fullfile(out_dir,sprintf('comparison.fig')));
saveas(h2,fullfile(out_dir,sprintf('comparison.png')));

h1=figure(4);
plot(freq,10*log10(S(:,7)));
ylabel('Power/Freq(dB/Hz)');
xlabel('Frequency (Hz)')
title('Using mtspectrumc');




