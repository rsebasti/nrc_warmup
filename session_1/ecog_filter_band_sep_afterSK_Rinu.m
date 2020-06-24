%% ECog band seperation
% Details of Data:
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Session : 1
%Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
%Assumption AC Freq: 60 Hz
% Directory -- creating directory if one does not exist.
base_dir=pwd;
echo off;
out_dir = fullfile(base_dir,'figs_ecog');
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end
warning('off','all')

% Data initializations and pre allocations
ecog_struct=load("bp_mot_t_h.mat");
ecog_data=ecog_struct.data;
clear ecog_struct;
[time_pts,channel_num]=size(ecog_data);
%Convert to microvolt
ecog_data= 0.0298*ecog_data;
%time vector for time in seconds
fs=1000;
t=(1:time_pts)/fs;
Fs = 1000;  % Sampling Frequency

% parameters for notch filter
Fnotch = 60;  % Notch Frequency
Q      = 20;  % Q-factor  % Bandwidth Attenuation
BW = Fnotch/Q;
[b, a] = iirnotch(Fnotch/(1000/2), BW/(1000/2), 1); % Calculate the coefficients using the IIRNOTCHPEAK function.

sep_freq=zeros(3,time_pts,channel_num);
yl = [min(min(ecog_data(:,1:channel_num))) max(max(ecog_data(:,1:channel_num)))];
filt_ecog_data=zeros([time_pts,channel_num]);
alpha_freq_band=zeros([time_pts,channel_num]);
beta_freq_band=zeros([time_pts,channel_num]);
gamma_freq_band=zeros([time_pts,channel_num]);

% Designing Filter
%1. First filter out 60 Hz from the signal for each channel
Hd = dsp.IIRFilter('Structure', 'Direct form II','Numerator', b,'Denominator', a);
%1) Alpha : 8 to 15 Hz; Bandpass at pass band freq - [8,15]
h1  = fdesign.bandpass('N,F3dB1,F3dB2', 6, 8, 15, 1000);
Hd1 = design(h1, 'butter');
%2) Beta : 16 to 31 Hz; Bandpass at pass band freq - [15,31]
h2  = fdesign.bandpass('N,F3dB1,F3dB2', 6, 16, 31, 1000);
Hd2 = design(h2, 'butter');
%3) Gamma : 32 Hz to 100 Hz, High Pass filter at Cut-off freq= 32 Hz
h3  = fdesign.bandpass('N,F3dB', 6, 32,100,1000);
Hd3 = design(h3, 'butter');
% Filtering Signal ( removing 60 Hz, segregating into alpha, beta and gamma
% bands)
for kk=1:channel_num
    % notching at 60 Hz
    filt_ecog_data(:,kk)=filtfilt(Hd.Numerator,Hd.Denominator,ecog_data(:,kk));
    % filtering alpha
    alpha_freq_band(:,kk)=filtfilt(Hd1.sosMatrix,Hd1.ScaleValues,filt_ecog_data(:,kk));
    % filtering beta
    beta_freq_band(:,kk)=filtfilt(Hd2.sosMatrix,Hd2.ScaleValues,filt_ecog_data(:,kk));
    % filtering gamma
    gamma_freq_band(:,kk)=filtfilt(Hd3.sosMatrix,Hd3.ScaleValues,filt_ecog_data(:,kk));
end
for kk=1:channel_num
    sep_freq(1,:,kk)=alpha_freq_band(:,kk);
    sep_freq(2,:,kk)=beta_freq_band(:,kk);
    sep_freq(3,:,kk)=gamma_freq_band(:,kk);
end
clear alpha_freq_band;
clear beta_freq_band;
clear gamma_freq_band;

y2 = [min(min(min(sep_freq))) max(max(max((sep_freq))))];

%Plotting
% 1. Plotting Sample Unflitered ecog signals
h0=figure;
% Creates a matrix of handles to your custom-spaced subplots.
rows = 10;
cols = 1;
% define spacing parameters
left_margin = 0.08; % left spacing on the figure
bottom_margin = 0.08; % bottom spacing on the figure
right_margin = 0.04; % right spacing on the figure
top_margin = 0.08; % top spacing on the figure
ax_spacing_horizontal = 0.02; % spacing between columns of axes
ax_spacing_vertical = 0.05; % spacing between rows of axes
ax_width = (1-left_margin-right_margin-(cols-1)*ax_spacing_vertical)/cols; % width of each axis
ax_height = (1-top_margin-bottom_margin-(rows-1)*ax_spacing_horizontal)/rows; % height of each axis
[r,c] = ndgrid(1:rows,1:cols);
ax = arrayfun(@(r,c)axes('Position',[left_margin+(c-1)*(ax_width+ax_spacing_vertical) bottom_margin+(r-1)*(ax_height+ax_spacing_horizontal) ax_width ax_height]),r(:),c(:),'UniformOutput',true);
for ii = rows:-1:1
    for jj = 1:cols
        plot(ax(ii,jj),t,ecog_data(:,rows-ii+1));
        ylab={['Ch ' num2str(rows-ii+1)],'Amp(\muV)'};
        ylim(ax(ii,jj),yl);
        xlim(ax(ii,jj),t([1 end]));
        set(ylabel(ax(ii,jj),ylab));
        if ii>1
            set(ax(ii,jj),'XTickLabel',[]);
        else 
            set(xlabel(ax(ii,jj),'Time (sec)'));
        end
    end
end
set(title("Sample Graph Amplitude vs time - Channels 1 to 10"));
saveas(h0,fullfile(out_dir,sprintf('Sample_AmpVsTime.fig')));
saveas(h0,fullfile(out_dir,sprintf('Sample_AmpVsTime.png')));

% 2. Plotting alpha, beta and gamma band
row1=3;col1=1;
[r1,c1] = ndgrid(1:row1,1:col1);
ax_width1 = (1-left_margin-right_margin-(col1-1)*ax_spacing_vertical)/col1; % width of each axis
ax_height1 = (1-top_margin-bottom_margin-(row1-1)*ax_spacing_horizontal)/row1; % height of each axis
y_labs={'Alpha Freq. Band ','Beta Freq. Band','Gamma Freq. Band'};
for kk=1:channel_num
    sub_dir = fullfile(out_dir,['Channel_' num2str(kk)]);
    if exist(sub_dir,'dir')~=7
        [tf,msg] = mkdir(sub_dir);
        assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
    end
    h1=figure;
    ax1 = arrayfun(@(r1,c1)axes('Position',[left_margin+(c1-1)*(ax_width1+ax_spacing_vertical) bottom_margin+(r1-1)*(ax_height1+ax_spacing_horizontal) ax_width1 ax_height1]),r1(:),c1(:),'UniformOutput',true);
    for mm=row1:-1:1
        for nn=1:col1
        plot(ax1(mm,nn),t,sep_freq(row1-mm+1,:,kk));
        ylim(ax1(mm,nn),y2);
        xlim(ax1(mm,nn),t([1 end]));
        ylab2={string(y_labs(row1-mm+1)),'Amp(\muV)'};
        set(ylabel(ax1(mm,nn),ylab2));              
        if mm>1
            set(ax1(mm,nn),'XTickLabel',[]);
        else 
            set(xlabel(ax1(mm,nn),'Time (sec)'));
        end
        end
    end
    hold on;
    title(['Channel:' num2str(kk)]);
    saveas(h1,fullfile(sub_dir,sprintf('Channel_%d.fig',kk)));
    saveas(h1,fullfile(sub_dir,sprintf('Channel_%d.png',kk)));
end
