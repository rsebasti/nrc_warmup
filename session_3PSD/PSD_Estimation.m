 %% ECog PSD Estimation : with and without windows!
% Ref for normalisation and multiplication with factor 2:
% Algorithm section of https://www.mathworks.com/examples/signal/mw/signal-ex20286697-periodogram
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
% AC Freq: 60 Hz

%Directory infos
base_dir=pwd;
data_dir='C:\Users\rseba\Documents\NRC\session_3PSD\mot_t_h\data';
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

%clipping 
Fnotch = 60; % AC Frequency to be notched.
freq= Fs*(0:(N/2))/N; %freq vector for PSD 
clip_freq=find(freq>1*Fnotch-5 & freq<1*Fnotch+5);
for aa=2:round((N/2)/Fnotch)
    temp = find(freq>(aa*Fnotch-5)& freq<(aa*Fnotch+5)); %indices corresponding to 60 Hz or its Harmonics +/- 5 Hz
    clip_freq=[clip_freq temp]; %frequency vector to be clipped.
end 

%without windows
%  zero padding & FFT calculation -> PSD estimation (FFT squared and normalized)
%clipping -> multiplying power of all freq with 2 except DC n Nyq. rate -> plotting

for kk=1:channel_num
    ecog_fft(:,kk)=fft(ecog_data(:,kk),N);
    psd=abs(ecog_fft(:,kk)).^2;
    pwr_temp =psd(1:N/2+1)/(time_pts*Fs); %throwing away neg freq. 
    pwr_temp(clip_freq)=NaN; %assigning NaN near 5Hz offset = clipping 
    Pwr_fft(:,kk)=[pwr_temp(1);2*pwr_temp(2:end-1);pwr_temp(end)]; % compensating for the power for neg freq. 
    [pxx(:,kk)] = periodogram(ecog_data(:,kk),[],N); % crosschecking with periodogram()
    pxx(clip_freq,kk)=NaN; %assigning NaN near 5Hz offset = clipping 
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
PSD_ham(clip_freq_1,:)=NaN;
PSD_han(clip_freq_1,:)=NaN;
pxx_welch_ham(clip_freq_1,:)=NaN;
pxx_welch_han(clip_freq_1,:)=NaN;


h=figure(1); % plotting non-windowed PSD
ax=define_position(2,1);
plot(ax(2,1),freq, 10*log10(Pwr_fft(:,7))+30);
ylabel(ax(2,1),{'PSD-FFT','Power/Freq(dB/Hz)'});
xticklabels(ax(2,1),[]);
ylim(ax(2,1),[-80 65]);
title('Channel 7');
plot(ax(1,1),freq,10*log10(pxx(:,7)));
ylim(ax(1,1),[-80 65]);
ylabel(ax(1,1),{'Periodogram','Power/Freq(dB/Hz)'});
xlabel(ax(1,1),'Frequency (Hz)')
saveas(h,fullfile(out_dir,sprintf('Channel_7_withoutWindows.fig')));
saveas(h,fullfile(out_dir,sprintf('Channel_7_withoutWindows.png')));

h1=figure(2); % plotting windowed PSD
ax1=define_position(4,1);
plot(ax1(4,1),f, 10*log10(PSD_ham(:,7)));
ylabel(ax1(4,1),{'PSD- Hamming','Power/Freq(dB/Hz)'});
xticklabels(ax1(4,1),[]);
ylim(ax1(4,1),[-40 65]);
title('Channel 7');
plot(ax1(3,1),f,10*log10(pxx_welch_ham(:,7)));
ylim(ax1(3,1),[-40 65]);
ylabel(ax1(3,1),{'pwelch- Hamming','Power/Freq(dB/Hz)'});
xticklabels(ax1(3,1),[]);
plot(ax1(2,1),f,10*log10(PSD_han(:,7)));
ylim(ax1(2,1),[-40 65]);
ylabel(ax1(2,1),{'PSD- Hanning','Power/Freq(dB/Hz)'});
xticklabels(ax1(2,1),[]);
plot(ax1(1,1),f,10*log10(pxx_welch_han(:,7)));
ylim(ax1(1,1),[-40 65]);
ylabel(ax1(1,1),{'pwelch- Hanning','Power/Freq(dB/Hz)'});
xlabel(ax1(1,1),'Frequency (Hz)')
saveas(h1,fullfile(out_dir,sprintf('Channel_7_Windows.fig')));
saveas(h1,fullfile(out_dir,sprintf('Channel_7_Windows.png')));

h2=figure(3); % plotting windowed PSD

ax2=define_position(3,1);
plot(ax2(3,1),freq, 30+10*log10(Pwr_fft(:,7)));
ylabel(ax2(3,1),{'PSD-FFT','Power/Freq(dB/Hz)'});
xticklabels(ax2(3,1),[]);
ylim(ax2(3,1),[-80 65]);
title('Channel 7');
plot(ax2(2,1),f,10*log10(PSD_ham(:,7)));
ylim(ax2(2,1),[-80 65]);
ylabel(ax2(2,1),{'PSD-Hamming','Power/Freq(dB/Hz)'});
xticklabels(ax2(2,1),[]);
plot(ax2(1,1),f,10*log10(PSD_han(:,7)));
ylim(ax2(1,1),[-80 65]);
ylabel(ax2(1,1),{'PSD- Hanning','Power/Freq(dB/Hz)'});
xlabel(ax2(1,1),'Frequency (Hz)')
saveas(h2,fullfile(out_dir,sprintf('comparison.fig')));
saveas(h2,fullfile(out_dir,sprintf('comparison.png')));





