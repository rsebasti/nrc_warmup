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
data_dir='C:\Users\rseba\Documents\NRC\session_4\mot_t_h\data';
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
[mat_psd_ham,mat_psd_han,PSD_ham,PSD_han,off_set,time]=specto(ecog_data,channel_num,time_pts,over_lap,Fs); %fn call to user defined fn to calculate PSD 
% crosschecking with pwelch()
for kk=1:channel_num
    pxx_welch_ham(:,kk) = pwelch(ecog_data(:,7),hamming(seg_length),off_set,nfft,Fs);
    pxx_welch_han(:,kk) = pwelch(ecog_data(:,7),hanning(seg_length),off_set,nfft,Fs);
end
PSD_ham(clip_freq_1,:)=NaN;
PSD_han(clip_freq_1,:)=NaN;
pxx_welch_ham(clip_freq_1,:)=NaN;
pxx_welch_han(clip_freq_1,:)=NaN;

h1=figure(2); % plotting windowed PSD
ax1=define_position(4,1);
plot(ax1(4,1),f,10*log10(PSD_ham(:,7)));
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

%Spectrogram
[S t f] = mtspecgramc(lfpdata, [.5 .25], params);
[spect,f,t]=spectrogram(ecog_data(:,7),1000,500,1024);
% h2=figure(4);
% imagesc(time,Fs*f/(2*3.14),log(abs(spect)));
% set(gca,'YDir', 'normal');
% ylabel('Frequency(Hz)');
% xlabel('Time(sec)');
% title('Using Spectrogram(): channel 7');
h3=figure(5);
imagesc(time,Fs*f/(2*3.14),log(abs(mat_psd_ham(:,:,7))));
set(gca,'YDir', 'normal');
ylabel('Frequency(Hz)');
xlabel('Time(sec)');
title('Using user defined fn: channel 7');
% h4=figure(6);
% imagesc(t,Fs*f/(2*3.14),(log(abs(spect))-log(abs(mat_psd_ham(:,:,7)))));
% set(gca,'YDir', 'normal');
% ylabel('Frequency(Hz)');
% title('Residual Spectrogram - diff(Spect() - userdefined()): channel 7');
% saveas(h2,fullfile(out_dir,sprintf('Sprectrogram.fig')));
% saveas(h2,fullfile(out_dir,sprintf('Sprectrogram.png')));
saveas(h3,fullfile(out_dir,sprintf('User_defined.fig')));
saveas(h3,fullfile(out_dir,sprintf('User_defined.png')));
% saveas(h4,fullfile(out_dir,sprintf('Residual.fig')));
% saveas(h4,fullfile(out_dir,sprintf('Residual.png')));



