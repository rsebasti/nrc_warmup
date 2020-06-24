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
% data_dir='C:\Users\rseba\Documents\NRC\session_4\mot_t_h\data';
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
N=2^nextpow2(time_pts); %Time point for FFT
ecog_fft=zeros([N,channel_num]);
Pwr_fft=zeros([N/2+1,channel_num]);
pxx=zeros([N/2+1,channel_num]);
over_lap=50; % 50% percentage overlap
seg_length=1000; %dividing signal into segments of 1000, essentially giving 1s segments
nfft=1024; 

%Spectrogram 3 ways: a) user_defined, b)using mtspecgramc() (Chronux), c)using
%spectrogram()
[mat_psd_ham,mat_psd_han,PSD_ham,PSD_han,off_set,time,freq]=specto(ecog_data,channel_num,time_pts,over_lap,Fs); %fn call to user defined fn to calculate PSD 
params.Fs=Fs;
[S,tr,fr] = mtspecgramc(ecog_data(:,7), [seg_length/Fs off_set/Fs], params);
[spect,f,t]=spectrogram(ecog_data(:,7),seg_length,off_set,nfft, Fs);
h2=figure(4);
imagesc(t,(Fs/2)*(f/pi),log(abs(spect)));
h21=colorbar;
ylabel(h21,'Power/Frequency (dB/Hz)')
set(gca,'YDir', 'normal');
ylabel('Frequency(Hz)');
xlabel('Time(s)');
title('Using spectrogram(): channel 7');
h3=figure(5);
imagesc(time,freq,log(abs(mat_psd_ham(:,:,7))));
h31=colorbar;
ylabel(h31,'Power/Frequency (dB/Hz)')
set(gca,'YDir', 'normal');
ylabel('Frequency(Hz)');
xlabel('Time(s)');
title('Using user defined fn: channel 7');
h4=figure(6);
imagesc(tr,fr,log(abs(S')))
h41=colorbar;
ylabel(h41,'Power/Frequency (dB/Hz)')
set(gca,'YDir', 'normal');
ylabel('Frequency(Hz)');
xlabel('Time(s)');
title('Using mtspecgramc(): channel 7');
saveas(h2,fullfile(out_dir,sprintf('using spectrogram.fig')));
saveas(h2,fullfile(out_dir,sprintf('using spectrogram.png')));
saveas(h3,fullfile(out_dir,sprintf('User_defined.fig')));
saveas(h3,fullfile(out_dir,sprintf('User_defined.png')));
saveas(h4,fullfile(out_dir,sprintf('using mtspecgramc.fig')));
saveas(h4,fullfile(out_dir,sprintf('using mtspecgramc.png')));



