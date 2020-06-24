%initializations
ecog=load('C:\Users\rseba\Documents\NRC\mot_t_h\mot_t_h\data\bp_mot_t_h.mat');
channel_num=size(ecog.data,2);
time=size(ecog.data,1);
figure(1);
for i=1:channel_num
     subplot(16,3,i)
     plot(ecog.data(:,i))
     title(['Channel: ' num2str(i)])
%      xlabel('Time(ms)');
%      ylabel('Amplitude(\muV)');
     hold on;
end