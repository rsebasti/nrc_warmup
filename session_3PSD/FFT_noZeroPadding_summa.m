%% ECog PSD Estimation : with and without windows!
% Ref for normalisation and multiplication with factor 2:
% Algorithm section of https://www.mathworks.com/examples/signal/mw/signal-ex20286697-periodogram
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
% AC Freq: 60 Hz

%Directory infos
base_dir=pwd;
data_dir='C:\Users\rseba\Documents\NRC\session_3\mot_t_h\data';
out_dir = fullfile(base_dir,'figs_PSD'); %creating directory for output if one does not exist
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end

%ecog data
ecog_struct=load(fullfile(data_dir,'bp_mot_t_h.mat'));
ecog_data=ecog_struct.data; 
clear ecog_struct;
ecog_data= 0.0298*ecog_data; %Convert to microvolt
[time_pts,channel_num]=size(ecog_data); % Size of ecog data matrix
Fs=1000;   % Sampling Frequency
t=(1:time_pts)/Fs; %time vector for time in seconds
N=2^nextpow2(time_pts); %Time point for FFT
freq= 1:((Fs/2)-1)/(N/2):Fs/2;%freq vector for PSD
ecog_fft=zeros([N,channel_num]);
Pwr_fft=zeros([N/2+1,channel_num]);

%clipping 
Fnotch = 60; % AC Frequency to be notched.
clip_freq=find(freq>1*Fnotch-5 & freq<1*Fnotch+5);
for aa=2:round((Fs/2)/Fnotch)
    temp = find(freq>aa*60-5 & freq<aa*60+5); %indices corresponding to 60 Hz or its Harmonics +/- 5 Hz
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
    Pwr_fft(:,kk)=[pwr_temp(1);2*pwr_temp(2:end)]; % compensating for the power for neg freq. 
end

%with windows
% Hamming window => element wise windowing-> zero padding & FFT calculation
% -> PSD estimation -> clipping -> plotting

ham_win=hamming(time_pts);
han_win=hanning(time_pts);
ecog_han=zeros([N,channel_num]);
ecog_ham=zeros([N,channel_num]);
ecog_han_fft=zeros([N,channel_num]);
ecog_ham_fft=zeros([N,channel_num]);
Pwr_ham=zeros([N/2+1,channel_num]);
Pwr_han=zeros([N/2+1,channel_num]);
pxx_ham=zeros([N/2+1,channel_num]);
pxx_han=zeros([N/2+1,channel_num]);
w_ham=zeros([N/2+1,channel_num]);
w_han=zeros([N/2+1,channel_num]);

for kk=1:channel_num
    ecog_han(1:time_pts,kk)=han_win.*ecog_data(1:time_pts,kk);
    ecog_ham(1:time_pts,kk)=ham_win.*ecog_data(1:time_pts,kk);
    ecog_ham_fft(:,kk)=fft(ecog_ham(:,kk),N);
    ecog_han_fft(:,kk)=fft(ecog_han(:,kk),N);
    psd_ham=abs(ecog_ham_fft(:,kk)).^2;
    psd_han=abs(ecog_han_fft(:,kk)).^2;  
    p_temp_ham =psd_ham(1:N/2+1)/(time_pts*Fs); %throwing away neg freq. 
    p_temp_han =psd_han(1:N/2+1)/(time_pts*Fs);
    p_temp_ham(clip_freq)=NaN;%assigning NaN for 5Hz offset around 60Hz and its harmonics.
    p_temp_han(clip_freq)=NaN;
    Pwr_ham(:,kk)=[p_temp_ham(1);2*p_temp_ham(2:end)];
    Pwr_han(:,kk)=[p_temp_han(1);2*p_temp_han(2:end)];
    %Using Periodogram()
    [pxx_ham(:,kk),w_ham(:,kk)] = periodogram(ecog_data(:,kk),hamming(time_pts));
    [pxx_han(:,kk),w_han(:,kk)] = periodogram(ecog_data(:,kk),hanning(time_pts));
    pxx_ham(clip_freq,kk)=NaN;
    pxx_han(clip_freq,kk)=NaN;
end

%plotting
rows=5;cols=1;
% define spacing parameters
left_margin = 0.08; % left spacing on the figure
bottom_margin = 0.08; % bottom spacing on the figure
right_margin = 0.04; % right spacing on the figure
top_margin = 0.08; % top spacing on the figure
ax_spacing_horizontal = 0.02; % spacing between columns of axes
ax_spacing_vertical = 0.01; % spacing between rows of axes
ax_width = (1-left_margin-right_margin-(cols-1)*ax_spacing_vertical)/cols; % width of each axis
ax_height = (1-top_margin-bottom_margin-(rows-1)*ax_spacing_horizontal)/rows; % height of each axis
[r,c] = ndgrid(1:rows,1:cols);


h=figure();
ax = arrayfun(@(r,c)axes('Position',[left_margin+(c-1)*(ax_width+ax_spacing_vertical) bottom_margin+(r-1)*(ax_height+ax_spacing_horizontal) ax_width ax_height]),r(:),c(:),'UniformOutput',true);
xticklabels(ax,[]);
yticklabels(ax,[]);
plot(ax(5,1),freq, 10*log10(Pwr_fft(:,7)));
ylabel(ax(5,1),{'PSD-FFT','Power/Freq(dB/Hz)'});
xticklabels(ax(5,1),[]);
ylim([-100 50]);
title('Channel 7');
plot(ax(4,1),freq,10*log10(Pwr_ham(:,7)));
ylabel(ax(4,1),{'PSD-Hamming','Power/Freq(dB/Hz)'});
xticklabels(ax(4,1),[]);
ylim([-100 50]);
plot(ax(3,1),freq,10*log10(Pwr_han(:,7)));
ylabel(ax(3,1),{'PSD-Hanning','Power/Freq(dB/Hz)'});
xticklabels(ax(3,1),[]);
ylim([-100 50]);
plot(ax(2,1),w_ham(:,7)*500/pi,10*log10(pxx_ham(:,7)));
ylabel(ax(2,1),{'PG-Hanning','Power/Freq(dB/Hz)'});
set(ax(2,1),'XTickLabel',[])
ylim([-100 50]);
plot(ax(1,1),w_han(:,7)*500/pi,10*log10(pxx_han(:,7)));
ylabel(ax(1,1),{'PG-Hamming','Power/Freq(dB/Hz)'});
xlabel(ax(1,1),'Frequency (Hz)')
ylim([-100 50]);

saveas(h,fullfile(out_dir,sprintf('Channel_7.fig')));
saveas(h,fullfile(out_dir,sprintf('Channel_7.png')));



