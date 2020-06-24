%initializations
%base_dir='C:\Users\rseba\Documents\NRC\mot_t_h\mot_t_h';
base_dir=pwd;
%ecog=load('C:\Users\rseba\Documents\NRC\mot_t_h\mot_t_h\data\bp_mot_t_h.mat');
ecog=load(pwd '\bp_mot_t_h.mat');
channel_num=size(ecog.data,2);
time=size(ecog.data,1);
mkdir figs
%plotting
for i=1:channel_num
     h=figure;
     plot(ecog.data(:,i))
     title(['Channel: ' num2str(i)])
     xlabel('Time(ms)');
     ylabel('Amplitude(\muV)');
     filename=strcat(base_dir,"\figs\Channel"+num2str(i)+".fig");
     savefig(h,filename);
end