%% Spectrogram
% Ref for normalisation and multiplication with factor 2:
% Algorithm section of https://www.mathworks.com/examples/signal/mw/signal-ex20286697-periodogram
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
% AC Freq: 60 Hz
close all;
%Directory infos
base_dir=pwd;
data_dir=pwd;
out_dir = fullfile(base_dir,'figs_PSD'); %creating directory for output if one does not exist
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end

%ecog data and initializations
ecog_struct=load(fullfile(data_dir,'bp_mot_t_h.mat'));
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

%with windows
% Hamming/ Hanning window => element wise windowing-> zero padding & FFT calculation
% -> PSD estimation -> clipping -> plotting
over_lap=50; % 50% percentage overlap
seg_length=1000; %dividing signal into segments of 1000, essentially giving 1s segments
nfft=1024; % 1024 point DFT
f=Fs*(0:(nfft/2))/nfft; %Freq vector for Windowed PSD calculation
% clipping
Fnotch = 60; % AC Frequency to be notched.
clip_freq_1=find(f>1*Fnotch-5 & f<1*Fnotch+5);
for aa=2:round((nfft/2)/Fnotch)
    temp = find(f>(aa*60-5)& f<(aa*60+5)); %indices corresponding to 60 Hz or its Harmonics +/- 5 Hz
    clip_freq_1=[clip_freq_1 temp]; %frequency vector to be clipped.
end 
pxx_welch_ham = zeros([nfft/2+1,channel_num]);
pxx_welch_han = zeros([nfft/2+1,channel_num]);
[mat_psd_ham,~,PSD_ham,~,off_set,time]=specto(ecog_data,channel_num,time_pts,over_lap,Fs); %fn call to user defined fn to calculate PSD 
% crosschecking with pwelch()
for kk=1:channel_num
    pxx_welch_ham(:,kk) = pwelch(ecog_data(:,7),hamming(seg_length),off_set,nfft,Fs);
end
PSD_ham(clip_freq_1,:)=NaN;
pxx_welch_ham(clip_freq_1,:)=NaN;


params.Fs=Fs;
params.tapers =[3 5];
M=1000; %1 second window
overlap_per=50;
N=2^nextpow2(M); 
off_set=overlap_per*M/100;
signal=ecog_data;
sum_Pxx_hamming_chronux=zeros(N/2+1,1);
k = fix((time_pts-off_set)/(M-off_set));
for kk=1:channel_num
    for ll=0:k-1
        signal_hamming=signal(1+ll*off_set:M+ll*off_set,kk).*hamming(M);
        [Spect,fr] = mtspectrumc(signal_hamming,params);
        sum_Pxx_hamming_chronux=sum_Pxx_hamming_chronux+Spect;
    end
    PSD_ham_chronux(:,kk)=sum_Pxx_hamming_chronux/k;
%     PSD_ham_chronux(2:end-1,kk) = PSD_ham(2:end-1,kk)*2;
%     PSD_han_chronux(:,kk)=sum_Pxx_hamming/(Fs*sum_hanning); 
%     PSD_han_chronux(2:end-1,kk) = PSD_han(2:end-1,kk)*2;
end
PSD_ham_chronux(clip_freq_1,:)=NaN;

h1=figure(1); % plotting windowed PSD
ax1=define_position(3,1);
plot(ax1(3,1),f,10*log10(PSD_ham(:,7)));
ylabel(ax1(3,1),{'PSD- Hamming','Power/Freq(dB/Hz)'});
xticklabels(ax1(3,1),[]);
plot(ax1(2,1),f,10*log10(pxx_welch_ham(:,7)));
ylabel(ax1(2,1),{'pwelch-Hamming','Power/Freq(dB/Hz)'});
xticklabels(ax1(2,1),[]);
plot(ax1(1,1),fr,10*log10(PSD_ham_chronux(:,7)));
ylabel(ax1(1,1),{'Chronux-Hamming','Power/Freq(dB/Hz)'});
xlabel(ax1(1,1),'Frequency (Hz)')
title('Channel 7 : User Defined, pwelch(), mtspectrumc()');
saveas(h1,fullfile(out_dir,sprintf('Channel_7_Windows.fig')));
saveas(h1,fullfile(out_dir,sprintf('Channel_7_Windows.png')));






