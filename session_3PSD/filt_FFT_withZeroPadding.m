%% ECog band seperation
% Details of Data:
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Session : 1
%Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
%Assumption AC Freq: 60 Hz
% Directory -- creating directory if one does not exist.
base_dir=pwd;
%creating directory for output
out_dir = fullfile(base_dir,'figs_PSD');
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end
data_dir='C:\Users\rseba\Documents\NRC\session_3\mot_t_h\data';
ecog_struct=load(fullfile(data_dir,'bp_mot_t_h.mat'));
ecog_data=ecog_struct.data;
clear ecog_struct;

[time_pts,channel_num]=size(ecog_data);
N=2^nextpow2(time_pts); %Time point for FFT
ecog_data= 0.0298*ecog_data;%Convert to microvolt

% ecog_data=[ecog_data;zeros(2^N-size(ecog_data,1),size(ecog_data,2))]; %zero padding

Fs=1000;  % Sampling Frequency
t=(1:time_pts)/Fs; %time vector for time in seconds
freq=0:Fs/time_pts:Fs/2; %freq vector for PSD
Fnotch = 60;  q = 20; bw = (Fnotch/(Fs/2))/q;
[b,a] = iircomb(round(Fs/Fnotch),bw,'notch');
clip_freq=find(freq>1*60-5 & freq<1*60+5);

for aa=2:round((Fs/2)/Fnotch)
    temp = find(freq>aa*60-5 & freq<aa*60+5);
    clip_freq=[clip_freq temp];
end
% [b, a] = iirnotch(Fnotch/(1000/2), BW/(1000/2), 1); % Calculate the coefficients using the IIRNOTCHPEAK function.
filt_ecog_data=zeros([time_pts,channel_num]);
filt=zeros(time_pts,1);
fft_ecog_data=zeros([N,channel_num]);
% fft_ecog_data=zeros([time_pts,channel_num]);
hann_ecog=zeros([time_pts,channel_num]);
hamm_ecog=zeros([time_pts,channel_num]);
% P=zeros([size(freq,2),channel_num]);
P=zeros([time_pts/2+1,channel_num]);
Hd = dsp.IIRFilter('Structure', 'Direct form II','Numerator', b,'Denominator', a);
for kk=1:channel_num
    % notching at 60 Hz
    filt_ecog_data(:,kk)=filtfilt(Hd.Numerator,Hd.Denominator,ecog_data(:,kk));
    fft_ecog_data(:,kk)=fft(filt_ecog_data(:,kk),N);
    psd=abs(fft_ecog_data(:,kk)).^2;
    p_temp =(time_pts*Fs)*psd(1:time_pts/2+1); %throwing away neg freq. 
    p_temp(clip_freq)=NaN;%assigning NaN near 5Hz offset
    P(:,kk)=[p_temp(1);2*p_temp(2:end-1);p_temp(end)];
end
h1 = figure('Position',[50 100 800 600]);
plot(freq, 10*log10(P(:,7)));
title('PSD Using notch filtered FFT with zeroPadding : Channel 7');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
saveas(h1,fullfile(out_dir,sprintf('PSD_Using_notch_filtered_FFT with_zeroPadding_Channel7.png')))
