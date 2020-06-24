%% EEG band seperation
% Details of Data:
% EEG Data: BCI Competition IV Dataset 1 ( 1000 Hz, calib), http://bbci.de/competition/iv/download/
% Reseach Group: Berlin Institute of Technology (Machine Learning Laboratory) and Fraunhofer FIRST (Intelligent Data Analysis Group) (Klaus-Robert Müller, Benjamin Blankertz, Carmen Vidaurre, Guido Nolte), and Campus Benjamin Franklin of the Charité - University Medicine Berlin, Department of Neurology, Neurophysics Group (Gabriel Curio) 
% No of Channels: 59
% Sampling Frequency : 1000 Hz
% Session : 1
% to convert to microV -> cnt= 0.1*double(cnt) where cnt is data -- ref: http://bbci.de/competition/iv/desc_1.html
%Signals were band-pass filtered between 0.05 and 200 Hz and then digitized at 1000 Hz with 16 bit (0.1 uV)
%Assumption AC Freq: 50 Hz (Germany)
clear all;
close all;

% Directory -- creating directory if one does not exist.
base_dir=pwd;
out_dir = fullfile(base_dir,'figs_assignment2');
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end

% Data initializations
eeg_struct=load("BCICIV_calib_ds1a_1000Hz.mat");
eeg_data=eeg_struct.cnt;
[time_pts,channel_num]=size(eeg_data);
%Convert to microvolt
eeg_data= 0.1*double(eeg_data);
yl=[min(eeg_data(:)) max(eeg_data(:))];
%time vector for time in seconds
fs=1000;
t=(1:time_pts)/fs;

%Filters used are Butterworth ( no ripples in passband ) and fifltfilt function to eliminate group delay 
% Finding filter coefficients for 6th order butterworth for 50 Hz band pass filter and for 
% delta, theta, alpha, beta, gamma bands for each channel

%1. First filter out Freq from 0-50Hz
[b_50,a_50] = butter(6,50/(fs/2));
%a) Delta : 0 to 4; Hz Low pass at Cut off @ 4Hz
[b_delta,a_delta] = butter(6,4/(fs/2));
%b) Theta : 4 to 8 Hz; Bandpass at pass band freq - [4,8]
[b_theta,a_theta] = butter(6,[4 8]/(fs/2));
%c) Alpha : 8 to 13 Hz; Bandpass at pass band freq - [8,13]
[b_alpha,a_alpha] = butter(6,[8 13]/(fs/2));
%d) Beta : 13 to 30 Hz; Bandpass at pass band freq - [13,30]
[b_beta,a_beta] = butter(6,[13 30]/(fs/2));
%d) Gamma : 30 Hz and above (until 50 Hz), Cut-off freq @ 30 Hz
[b_gamma,a_gamma] = butter(6,30/(fs/2),'high');

for ii=1:channel_num
    filt_eeg(:,ii) = filter(b_50,a_50,eeg_data(:,ii));
    %Filtfilt : for alpha, beta, gamma, delta, theta bands
    eeg_delta(:,ii) = filtfilt(b_delta,a_delta,filt_eeg(:,ii));
    eeg_theta(:,ii) = filtfilt(b_theta,a_theta,filt_eeg(:,ii));
    eeg_alpha(:,ii) = filtfilt(b_alpha,a_alpha,filt_eeg(:,ii));
    eeg_beta(:,ii) = filtfilt(b_beta,a_beta,filt_eeg(:,ii));
    eeg_gamma(:,ii) = filtfilt(b_gamma,a_gamma,filt_eeg(:,ii));
end

%plotting
for ii=1:channel_num
    % Creating a subdirectory for each channel and plots related to them.
    sub_dir = fullfile(out_dir,['Channel_' num2str(ii)]);
    if exist(sub_dir,'dir')~=7
        [tf,msg] = mkdir(sub_dir);
        assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
    end
    %Plotting and Saving Amplitude vs Time plot for each channel
    h1=figure(1);
    plot(t,eeg_data(:,ii))
    title(['Channel: ' num2str(ii)])
    xlabel('Time (sec)');
    ylabel('Amplitude (\muV)');
    xlim(t([1 end]));
    ylim(yl);
    saveas(h1,fullfile(sub_dir,sprintf('Channel%02d.fig',ii)));
    saveas(h1,fullfile(sub_dir,sprintf('Channel%02d.png',ii)));
    h2=figure(2);
    subplot(1,3,1)
    plot(t,eeg_delta(:,ii))
    title(['Delta - Channel: ' num2str(ii)])
    xlabel('Time (sec)');
    ylabel('Amplitude (\muV)');
    xlim(t([1 end]));
    ylim([min(eeg_delta(:,ii)) max(eeg_delta(:,ii))]);
    hold on;
%         plot(t,eeg_theta(:,ii))
%         title(['Theta - Channel: ' num2str(ii)])
%         xlabel('Time (sec)');
%         ylabel('Amplitude (\muV)');
%         xlim(t([1 end]));
%         ylim([min(eeg_theta(:,ii)) max(eeg_theta(:,ii))]);
%         hold on;
%         plot(t,eeg_alpha(:,ii))
%         title(['Alpha - Channel: ' num2str(ii)])
%         xlabel('Time (sec)');
%         ylabel('Amplitude (\muV)');
%         xlim(t([1 end]));
%         ylim([min(eeg_alpha(:,ii)) max(eeg_alpha(:,ii))]);
%         hold on;
    subplot(1,3,2)
    plot(t,eeg_beta(:,ii))
    title(['Beta - Channel: ' num2str(ii)])
    xlabel('Time (sec)');
    ylabel('Amplitude (\muV)');
    xlim(t([1 end]));
    ylim([min(eeg_beta(:,ii)) max(eeg_beta(:,ii))]);
    hold on;
    subplot(1,3,3)
    plot(t,eeg_gamma(:,ii))
    title(['Gamma - Channel: ' num2str(ii)])
    xlabel('Time (sec)');
    ylabel('Amplitude (\muV)');
    xlim(t([1 end]));
    ylim([min(eeg_gamma(:,ii)) max(eeg_gamma(:,ii))]);
    hold on;
    saveas(h2,fullfile(sub_dir,sprintf('Bands_Channel%02d.fig',ii)));
    saveas(h2,fullfile(sub_dir,sprintf('Bands_Channel%02d.png',ii)));
    close all;
end





