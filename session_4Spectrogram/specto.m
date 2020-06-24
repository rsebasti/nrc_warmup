function [mat_psd_ham,mat_psd_han,PSD_ham,PSD_han,off_set,time]=specto(signal,channel_num,time_pts,overlap_per,Fs)
M=1000; %1 second window
N=2^nextpow2(M); 
off_set=overlap_per*M/100;
% t_pts=size(signal,2);
f_pts=N/2+1;
k = fix((time_pts-off_set)/(M-off_set));
coloffsets = (0:(k-1))*(M-off_set);
time = (coloffsets+(M/2)')/Fs;
%L= fix((t_pts-off_set)/(length(window-noverlap)));
mat_psd_han=zeros(f_pts,k,channel_num);
mat_psd_ham=zeros(f_pts,k,channel_num);
sum_hamming=sum(hamming(M).^2);
sum_hanning=sum(hanning(M).^2);
PSD_ham=zeros(N/2+1,channel_num);
PSD_han=zeros(N/2+1,channel_num);
sum_Pxx_hamming=zeros(N/2+1,1);
sum_Pxx_hanning=zeros(N/2+1,1);
for kk=1:channel_num
    for ll=0:k-1
        signal_hamming=signal(1+ll*off_set:M+ll*off_set,kk).*hamming(M);
        signal_hanning=signal(1+ll*off_set:M+ll*off_set,kk).*hanning(M);
        fft_hamming=fft(signal_hamming,N);
        fft_hanning=fft(signal_hanning,N);
        psd_hamming=abs(fft_hamming).^2;
        psd_hanning=abs(fft_hanning).^2;
        mat_psd_ham(:,ll+1,kk)=psd_hamming(1:N/2+1)/(Fs*sum_hamming);
        mat_psd_han(:,ll+1,kk)=psd_hanning(1:N/2+1)/(Fs*sum_hanning); %throwing away neg freq.
        pwr_temp_hamming=psd_hamming(1:N/2+1);
        pwr_temp_hanning=psd_hanning(1:N/2+1); %throwing away neg freq.
        sum_Pxx_hamming=sum_Pxx_hamming+pwr_temp_hamming;
        sum_Pxx_hanning=sum_Pxx_hanning+pwr_temp_hanning;
    end
    PSD_ham(:,kk)=sum_Pxx_hamming/(Fs*sum_hamming);
    PSD_ham(2:end-1,kk) = PSD_ham(2:end-1,kk)*2;
    PSD_han(:,kk)=sum_Pxx_hamming/(Fs*sum_hanning); 
    PSD_han(2:end-1,kk) = PSD_han(2:end-1,kk)*2;
end

    

    
    
        


