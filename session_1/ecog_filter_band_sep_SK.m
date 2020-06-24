%% EEG band seperation
% Details of Data:
% No of Channels: 47
% Sampling Frequency : 1000 Hz
% Session : 1
%Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
%Assumption AC Freq: 60 Hz
clear all;
close all;
% clf; % "clf" (clear figure) opens up an unused figure you don't need; you
% can use "clc" (clear screen) if you want to reset your command window
% display
clc;
% as a side note, I generally don't use clear all/close all/clc any more.
% if there's really a situation where I want/need it, I'll put that code
% in a function (which starts with its own blank workspace).

% Directory -- creating directory if one does not exist.
base_dir=pwd;
out_dir = fullfile(base_dir,'figs_ecog');
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end
warning('off','all')

% Data initializations and pre allocations
data_dir = 'C:\Users\rseba\Documents\NRC\session1'; % specify location of data file
ecog_struct=load(fullfile(data_dir,'bp_mot_t_h.mat'));
ecog_data=ecog_struct.data;
clear ecog_struct; % clear unused struct to save memory
[time_pts,channel_num]=size(ecog_data);
%Convert to microvolt
ecog_data= 0.0298*ecog_data;
%time vector for time in seconds
fs=1000;
t=(1:time_pts)/fs;
yl = [min(min(ecog_data(:,1:channel_num))) max(max(ecog_data(:,1:channel_num)))];
notch_ecog_data=zeros([time_pts,channel_num]);
% ecog_data_temp=zeros([time_pts,channel_num]); % since you're not going to
% use ecog_data_temp after processing any given channel, it doesn't make
% sense to allocate and fill up the space to copy all the data across all
% channels
%alpha_band_temp=zeros([time_pts,channel_num]); % similarly for all the
%"temps" below - saw you were applying filter, then filtfilt - this isn't
%necessary (unless I misunderstand what you're trying to do). just need the
%filtfilt.
alpha_freq_band=zeros([time_pts,channel_num]);
%beta_band_temp=zeros([time_pts,channel_num]); 
beta_freq_band=zeros([time_pts,channel_num]);
%gamma_band_temp=zeros([time_pts,channel_num]); 
gamma_freq_band=zeros([time_pts,channel_num]);

% Designing Filter
%1. First filter out 60 Hz from the signal for each channel
h  = fdesign.bandstop('N,F3dB1,F3dB2', 6, 58, 62, 1000);
Hd = design(h, 'butter');
%1) Alpha : 8 to 15 Hz; Bandpass at pass band freq - [8,15]
h1  = fdesign.bandpass('N,F3dB1,F3dB2', 6, 8, 15, 1000);
Hd1 = design(h1, 'butter');
%2) Beta : 16 to 31 Hz; Bandpass at pass band freq - [16,32]
h2  = fdesign.bandpass('N,F3dB1,F3dB2', 6, 16, 32, 1000);
Hd2 = design(h2, 'butter');
%3) Gamma : 32 Hz to 100 Hz, Bandpass at pass band freq - [32,100] - don't
%want to keep all the high frequency noisy stuff here if you don't need it
h3  = fdesign.bandpass('N,F3dB1,F3dB2', 6, 32, 100, 1000);
Hd3 = design(h3, 'butter');
% Filtering Signal ( removing 60 Hz, segregating into alpha, beta and gamma
% bands)
for kk=1:channel_num
    
    %ecog_data_temp=filter(Hd,ecog_data(:,kk));
    
    % not certain why filter, then filtfilt, with the same Hd object;
    % probably only need one (the filtfilt)
    notch_ecog_data(:,kk)=filtfilt(Hd.sosMatrix,Hd.ScaleValues,ecog_data(:,kk));
    
    % filtering alpha
    %alpha_band_temp(:,kk)=filter(Hd1,filt_ecog_data(:,kk));
    alpha_freq_band(:,kk)=filtfilt(Hd1.sosMatrix,Hd1.ScaleValues,ecog_data(:,kk));
    
    % filtering beta
    %beta_band_temp(:,kk)=filter(Hd2,filt_ecog_data(:,kk));
    beta_freq_band(:,kk)=filtfilt(Hd2.sosMatrix,Hd2.ScaleValues,ecog_data(:,kk));
    
    % filtering gamma
    %gamma_band_temp(:,kk)=filter(Hd3,filt_ecog_data(:,kk));
    gamma_freq_band(:,kk)=filtfilt(Hd3.sosMatrix,Hd3.ScaleValues,ecog_data(:,kk));
end
% this structure: min([data(:) data(:) data(:)]) uses a lot of memory
% because you have to create a big matrix that contains all of your data at
% once. If you instead do [min(data(:)) min(data(:)) min(data(:))] you can
% conserve a lot of memory and probably some processing time.
y2 = [min([min(alpha_freq_band(:)) min(beta_freq_band(:)) min(gamma_freq_band(:))]) ... 
    max([max(alpha_freq_band(:)) max(beta_freq_band(:)) max(gamma_freq_band(:))])];

%Plotting
% 1. Plotting Sample Unflitered ecog signals
h0=figure; % don't need the (1) unless you want to refer to the figure like that later
% Creates a matrix of handles to your custom-spaced subplots.
rows = 10;
cols = 1;
space = 0.09;
% ax = getCustomAxesPos(rows,cols,space); % I'm going to suggest that you
% learn to do this from scratch on your own. ;) I've got some sample code
% below. You'll need to adjust the parameters, and make sure you go through
% it to understand what's going on. As an exercise, see if you can change
% the vertical order of the plots (top axes- channel 1, bottom axes-
% channel 2). Then, implement the same thing for the second plot (with the
% frequency bands).

% define spacing parameters
left_margin = 0.08; % left spacing on the figure
bottom_margin = 0.08; % bottom spacing on the figure
right_margin = 0.04; % right spacing on the figure
top_margin = 0.13; % top spacing on the figure
ax_spacing_horizontal = 0.02; % spacing between columns of axes
ax_spacing_vertical = 0.05; % spacing between rows of axes
ax_width = (1-left_margin-right_margin-(cols-1)*ax_spacing_vertical)/cols; % width of each axis
ax_height = (1-top_margin-bottom_margin-(rows-1)*ax_spacing_horizontal)/rows; % height of each axis

% potentially 2d array of axes objects - in this case you can specify the
% element in the "last" position (rightmost,bottommost) and that will
% pre-allocate the whole array to an array of axes objects.
[r,c] = ndgrid(1:rows,1:cols);
ax = arrayfun(@(r,c)axes('Position',[left_margin+(c-1)*(ax_width+ax_spacing_vertical) bottom_margin+(r-1)*(ax_height+ax_spacing_horizontal) ax_width ax_height]),r(:),c(:),'UniformOutput',true);
for ii = 1:rows
    for jj = 1:cols
        % plot data into the axes
        plot(ax(ii,jj),t,ecog_data(:,ii));
        % y-label for this axes
        ylab={['Channel ' num2str(ii)],'Amplitude(\muV)'};
        ylabel(ax(ii,jj),ylab);
        % x/y lims
        ylim(ax(ii,jj),yl);
        xlim(ax(ii,jj),t([1 end]));
        % font size
        % set(ax(ii,jj),'fontsize',3); % never put small fonts on figures! :)
        % set(ylabel(ax(ii,jj),ylab), 'fontsize', 4);
        if ii>1
            set(ax(ii,jj), 'XTick',[], 'XTickLabel',[]); % you need to set "XTick" to empty to get rid of the labels
        else
            xlabel(ax(ii,jj),'Time (sec)');
        end
    % hold on; % not sure why this hold on...
    end
end

% Only ever make font sizes bigger :)
title(ax(end,end),'Sample Graph Amplitude vs time (Channels 1-10)');
%set(title("Sample Graph Amplitude vs time - Channels 1 to 10"),'fontsize',9);

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



