%initializations
%base_dir='C:\Users\rseba\Documents\NRC\mot_t_h\mot_t_h';
base_dir=pwd;
src_file='bp_mot_t_h.mat';
%ecog=load('C:\Users\rseba\Documents\NRC\mot_t_h\mot_t_h\data\bp_mot_t_h.mat');
ecog=load(fullfile(pwd,src_file));
channel_num=size(ecog.data,2);
time=size(ecog.data,1);
mkdir figs_R
%plotting
for i=1:channel_num
     h=figure;
     plot(ecog.data(:,i))
     title(['Channel: ' num2str(i)])
     xlabel('Time(ms)');
     ylabel('Amplitude(\muV)');
     filename=strcat(base_dir,"\figs_R\Channel"+num2str(i)+".fig");
     savefig(h,filename);
end