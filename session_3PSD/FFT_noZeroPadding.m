%% ECog band seperation
% Details of Data:
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Session : 3
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
fft_ecog_data=zeros([time_pts,channel_num]);

P=zeros([time_pts/2+1,channel_num]);
Hd = dsp.IIRFilter('Structure', 'Direct form II','Numerator', b,'Denominator', a);
for kk=1:channel_num
    % notching at 60 Hz
%     filt_ecog_data(:,kk)=filtfilt(Hd.Numerator,Hd.Denominator,ecog_data(:,kk));
    fft_ecog_data(:,kk)=fft(ecog_data(:,kk));
    psd=abs(fft_ecog_data(:,kk)).^2;
    p_temp =(time_pts*Fs)*psd(1:time_pts/2+1); %throwing away neg freq. 
    p_temp(clip_freq)=NaN;%assigning NaN near 5Hz offset
    P(:,kk)=[p_temp(1);2*p_temp(2:end-1);p_temp(end)];
end
% yl = [min(min(filt_ecog_data(:,1:channel_num))) max(max(filt_ecog_data(:,1:channel_num)))];
h1 = figure('Position',[50 100 800 600]); %
% modified=[freq; 10*log10(P(:,7))];
plot(freq, 10*log10(P(:,7)));
title('PSD Using FFT without zeroPadding : Channel 7');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
saveas(h1,fullfile(out_dir,sprintf('PSD_FFT without_zeroPadding_Channel7.png')))

%with windows
% Hamming window => element wise windowing-> zero padding-> FFT calculation
% -> PSD estimation -> clipping -> plotting

ham_win=hamming(time_pts);
han_win=hanning(time_pts);

ecog_han=zeros([N,channel_num]);
ecog_ham=zeros([N,channel_num]);
ecog_han_fft=zeros([N,channel_num]);
ecog_ham_fft=zeros([N,channel_num]);
P_ham=zeros([N/2+1,channel_num]);
P_han=zeros([N/2+1,channel_num]);
freq= 1:499/262144:500;

clip_freq=find(freq>1*60-5 & freq<1*60+5);
for aa=2:round((Fs/2)/Fnotch)
    temp = find(freq>aa*60-5 & freq<aa*60+5);
    clip_freq=[clip_freq temp];
end

for kk=1:channel_num
    ecog_han(1:time_pts,kk)=han_win.*ecog_data(1:time_pts,kk);
    ecog_ham(1:time_pts,kk)=ham_win.*ecog_data(1:time_pts,kk);
    ecog_ham_fft(:,kk)=fft(ecog_ham(:,kk),N);
    ecog_han_fft(:,kk)=fft(ecog_han(:,kk),N);
    psd_ham=abs(ecog_ham_fft(:,kk)).^2;
    psd_han=abs(ecog_han_fft(:,kk)).^2;
    
    p_temp_ham =(time_pts*Fs)*psd_ham(1:N/2+1); %throwing away neg freq. 
    p_temp_han =(time_pts*Fs)*psd_han(1:N/2+1);
    p_temp_ham(clip_freq)=NaN;%assigning NaN for 5Hz offset around 60Hz and its harmonics.
    p_temp_han(clip_freq)=NaN;
    P_ham(:,kk)=[p_temp_ham(1);2*p_temp_ham(2:end-1);p_temp_ham(end)];
    P_han(:,kk)=[p_temp_han(1);2*p_temp_han(2:end-1);p_temp_han(end)];
    
end


% figure();
% subplot(2,1,1)
% plot(ecog_ham(:,7));
% subplot(2,1,2)
% plot(ecog_han(:,7));

h2 = figure('Position',[50 100 800 600]); %
plot(freq,10*log10(P_ham(:,7)));
title('PSD Using FFT (Hamming window): Channel 7');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
saveas(h2,fullfile(out_dir,sprintf('PSD_FFT_Hamming_Ch7.png')))

h3 = figure('Position',[50 100 800 600]); %
plot(freq,10*log10(P_han(:,7)));
title('PSD Using FFT (Hanning window): Channel 7');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
saveas(h3,fullfile(out_dir,sprintf('PSD_FFT_Hanning_Ch7.png')))
% 
% %Cross checking
% % figure();
% % [pxx_hamming,w_hamming] = periodogram(ecog_data(:,7),hamming(length(ecog_data(:,7))));
% % plot(w_hamming*500/pi,10*log10(pxx_hamming));
% % 
% % figure();
% % [pxx_han,w_han] = periodogram(ecog_data(:,7),hanning(length(ecog_data(:,7))));
% % plot(w_han*500/pi,10*log10(pxx_han));
% 
% 
% 
