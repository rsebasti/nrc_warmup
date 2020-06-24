function [PSD_ham,PSD_han,off_set]=PSD_calc(signal,channel_num,time_pts,overlap_per,Fs)
k=ceil(time_pts/Fs);
M=1000; %1 second window
N=2^nextpow2(M); 
off_set=overlap_per*M/100;
sum_hamming=sum(hamming(1000).^2);
sum_hanning=sum(hanning(1000).^2);
PSD_ham=zeros(N/2+1,channel_num);
PSD_han=zeros(N/2+1,channel_num);
sum_Pxx_hamming=zeros(N/2+1,1);
sum_Pxx_hanning=zeros(N/2+1,1);
for kk=1:channel_num
    for ll=0:k-1
        signal_hamming=signal(1+ll*off_set:M+ll*off_set,kk).*hamming(1000);
        signal_hanning=signal(1+ll*off_set:M+ll*off_set,kk).*hanning(1000);
        fft_hamming=fft(signal_hamming,N);
        fft_hanning=fft(signal_hanning,N);
        psd_hamming=abs(fft_hamming).^2;
        psd_hanning=abs(fft_hanning).^2;
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

    
    
        


