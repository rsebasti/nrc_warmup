%initializations
%base_dir='C:\Users\rseba\Documents\NRC\mot_t_h\mot_t_h';
base_dir=pwd;
src_file='bp_mot_t_h.mat';
out_dir = fullfile(base_dir,'figs');
if exist(out_dir,'dir')~=7
    [tf,msg] = mkdir(out_dir);
    assert(tf,'Could not create output directory "%s": %s',out_dir,msg);
end

% load data
%ecog=load('C:\Users\rseba\Documents\NRC\mot_t_h\mot_t_h\data\bp_mot_t_h.mat');
ecog=load(fullfile(pwd,src_file));
channel_num=size(ecog.data,2);

% construct time vector
fs = 1000;
t=(1:size(ecog.data,1))/fs;

% scale back into microvolts
yscale = 0.0298;
ecog.data = ecog.data*yscale; % convert to microvolts
yl = [min(ecog.data(:)) max(ecog.data(:))];

%plotting
for ii=1:channel_num
     h=figure;
     plot(t,ecog.data(:,ii))
     title(['Channel: ' num2str(ii)])
     xlabel('Time (sec)');
     ylabel('Amplitude (\muV)');
     xlim(t([1 end]));
     ylim(yl);
     saveas(h,fullfile(out_dir,sprintf('Channel%02d.fig',ii)));
     saveas(h,fullfile(out_dir,sprintf('Channel%02d.png',ii)));
     %savefig(h,filename);
end