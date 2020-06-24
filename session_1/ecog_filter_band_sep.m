%% ECog band seperation
% Details of Data:
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Session : 1
%Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
%Assumption AC Freq: 60 Hz
clear all;
close all;
% Directory -- creating directory if one does not exist.
base_dir=pwd;
out_dir = fullfile(base_dir,'figs_ecog');
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end
warning('off','all')

% Data initializations and pre allocations
ecog_struct=load("bp_mot_t_h.mat");
ecog_data=ecog_struct.data;
[time_pts,channel_num]=size(ecog_data);
%Convert to microvolt
ecog_data= 0.0298*ecog_data;
%time vector for time in seconds
fs=1000;
t=(1:time_pts)/fs;
yl = [min(min(ecog_data(:,1:channel_num))) max(max(ecog_data(:,1:channel_num)))];
filt_ecog_data=zeros([time_pts,channel_num]);
ecog_data_temp=zeros([time_pts,channel_num]);
alpha_band_temp=zeros([time_pts,channel_num]);
alpha_freq_band=zeros([time_pts,channel_num]);
beta_band_temp=zeros([time_pts,channel_num]); 
beta_freq_band=zeros([time_pts,channel_num]);
gamma_band_temp=zeros([time_pts,channel_num]); 
gamma_freq_band=zeros([time_pts,channel_num]);

% Designing Filter
%1. First filter out 60 Hz from the signal for each channel
h  = fdesign.bandstop('N,F3dB1,F3dB2', 6, 58, 62, 1000);
Hd = design(h, 'butter');
%1) Alpha : 8 to 15 Hz; Bandpass at pass band freq - [8,15]
h1  = fdesign.bandpass('N,F3dB1,F3dB2', 6, 8, 15, 1000);
Hd1 = design(h1, 'butter');
%2) Beta : 16 to 31 Hz; Bandpass at pass band freq - [16,31]
h2  = fdesign.bandpass('N,F3dB1,F3dB2', 6, 16, 31, 1000);
Hd2 = design(h2, 'butter');
%3) Gamma : 32 Hz and above, High Pass filter at Cut-off freq= 32 Hz
h3  = fdesign.highpass('N,F3dB', 6, 32, 1000);
Hd3 = design(h3, 'butter');
% Filtering Signal ( removing 60 Hz, segregating into alpha, beta and gamma
% bands)
for kk=1:channel_num
    ecog_data_temp(:,kk)=filter(Hd,ecog_data(:,kk));
    filt_ecog_data(:,kk)=filtfilt(Hd.sosMatrix,Hd.ScaleValues,ecog_data_temp(:,kk));
    % filtering alpha
    alpha_band_temp(:,kk)=filter(Hd1,filt_ecog_data(:,kk));
    alpha_freq_band(:,kk)=filtfilt(Hd1.sosMatrix,Hd1.ScaleValues,alpha_band_temp(:,kk));
    % filtering beta
    beta_band_temp(:,kk)=filter(Hd2,filt_ecog_data(:,kk));
    beta_freq_band(:,kk)=filtfilt(Hd2.sosMatrix,Hd2.ScaleValues,beta_band_temp(:,kk));
    % filtering gamma
    gamma_band_temp(:,kk)=filter(Hd3,filt_ecog_data(:,kk));
    gamma_freq_band(:,kk)=filtfilt(Hd3.sosMatrix,Hd3.ScaleValues,gamma_band_temp(:,kk));
end
y2 = [(min(min([alpha_freq_band(:) beta_freq_band(:) gamma_freq_band(:)]))) ...
    max(max([alpha_freq_band(:) beta_freq_band(:) gamma_freq_band(:)]))];

%Plotting
% 1. Plotting Sample Unflitered ecog signals
h0=figure(1);
% Creates a matrix of handles to your custom-spaced subplots.
rows = 10;
cols = 1;
space = 0.09;
ax = getCustomAxesPos(rows,cols,space);
for ii = 1:rows
    for jj = 1:cols
        plot(ax(ii,jj),t,ecog_data(:,ii));
        ylab={['Channel ' num2str(ii)],'Amplitude(\muV)'};
        ylim(ax(ii,jj),yl);
        xlim(ax(ii,jj),t([1 end]));
        set(ax(ii,jj),'fontsize',3);
        set(ylabel(ax(ii,jj),ylab), 'fontsize', 4);
        if ii<rows
            set(ax(ii,jj), 'XTickLabel',[]);
        else
            set(xlabel(ax(ii,jj),'Time (sec)'), 'fontsize', 9);
        end
    end
end
set(title("Sample Graph Amplitude vs time - Channels 1 to 10"),'fontsize',9);
saveas(h0,fullfile(out_dir,sprintf('Sample_AmpVsTime.fig')));
saveas(h0,fullfile(out_dir,sprintf('Sample_AmpVsTime.png')));
% 2. Plotting alpha, beta and gamma band
for ll=1:channel_num
    sub_dir = fullfile(out_dir,['Channel_' num2str(ll)]);
    if exist(sub_dir,'dir')~=7
        [tf,msg] = mkdir(sub_dir);
        assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
    end
    clf;
    h2=figure(2);
    ax1=subplot(3,1,1);
    plot(t,alpha_freq_band(:,ll));
    ylabel({'Alpha Freq Band','Amplitude(\muV)'});
    hold on;
    ax2=subplot(3,1,2);
    plot(t,beta_freq_band(:,ll));
    ylabel({'Beta Freq Band','Amplitude(\muV)'});
    hold on;
    ax3=subplot(3,1,3);
    plot(t,gamma_freq_band(:,ll));
    ylabel({'Gamma Freq Band','Amplitude(\muV)'});
    hold on;
    currentFigure = gcf;
    title(currentFigure.Children(end), ['Channel:' num2str(ll)]);
    linkaxes([ax1,ax2,ax3],'xy');
    xlim(ax1,t([1 end]));
    ylim(ax1,y2);
    set([ax1,ax2], 'XTickLabel',[]);
    xlabel(ax3,'Time (sec)');
    saveas(h2,fullfile(sub_dir,sprintf('Channel_%d.fig',ll)));
    saveas(h2,fullfile(sub_dir,sprintf('Channel_%d.png',ll)));
end



