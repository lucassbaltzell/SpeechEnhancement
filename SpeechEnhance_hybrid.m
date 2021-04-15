function [y] = SpeechEnhance_hybrid(stim,fs,ear,Tfast,Tslow_exp,Fctf,Tslow_cmp,gain_bounds_dB)
%This function takes a stereo input and ouputs an enhanced signal.
%stim: should be stereo input
%fs: sample rate of input
%ear: determines whether enhancement is independent across ears (ear =
%0),or whether ehnahement is linked across ears using the left (ear = 1) or
%right (ear = 2) as reference
%Tfast: short-duration time frame (in seconds) over which gain is applied. 
%Becuase this is an FFT-based implementation, this value will be rounded up 
%to a duration corresponding to the nearest power of 2. Because we use a 50 
%percent overlap-and-add procedure (windowed), this refers to the "effective" 
%time frame, while the actual time frame has a durataion twice that of Tfast.
%Tslow_exp: long-duration time frame (in seconds) for normalized expansive 
%gain. Normalization refers to the fact that the maximum expansive gain
%that can be applied is 1. In other words, expansion has the effect of
%driving low-intensity time frames closer to zero.
%Fctf: Frequency cutoff for centroid. Determines boundary between voiceless
%"consonant" bins and voiced "vowel" bins
%Tslow_cmp: long-duration time frame (in seconds) for compressive gain.
%gain_bounds_dB = limits for compressive gain
%This compressive gain is similar to the EEQ method described in 
%Desloge et al (2017) JASA, with corresponding Patent: US20170208399A1

%[y] = SpeechEnhance_hybrid(stim,fs,ear,0.005,0.05,1500,0.5,[0 20])

%created by Lucas Baltzell 04/15/21

if nargin == 2
    ear = 0;
    Tfast = 0.004; %in seconds, and corresponding to 64-point FFT 
    Tslow_exp = 0.05; %s
    Fctf = 1500; %Hz
    Tslow_cmp = 0.5; %s
    gain_bounds_dB = [0 20];
elseif nargin == 3
    Tfast = 0.004;
    Tslow_exp = 0.05;
    Fctf = 1500;
    Tslow_cmp = 0.5;
    gain_bounds_dB = [0 20];
end
    
if fs ~= 16000
    stim = resample(stim,16000,fs);
    fs = 16000;
end

target_step = Tfast*fs;
NFFT = 2.^round(log2(2*target_step)); %NFFT is twice Tfast for 50% overlap/add
[XL,f] = spectrogram(cat(1,stim(:,1),zeros(NFFT/2,1)),boxcar(NFFT),NFFT/2,NFFT,fs);
XR = spectrogram(cat(1,stim(:,2),zeros(NFFT/2,1)),boxcar(NFFT),NFFT/2,NFFT,fs);
nt = size(XL,2);

nbin_exp = round(Tslow_exp/Tfast);
nbin_cmp = round(Tslow_cmp/Tfast);
%exp params%
tau = 0.5;
kv = 2;
kc = 0.1;

%cmp params%
gain_bounds_A = 10.^(gain_bounds_dB/20);

XT_enh_L = [];
XT_enh_R = [];
for n = 1:nt
    %%run expansion%%
    t0 = max(1,n-nbin_exp+1);
    xwinL = XL(:,t0:n);
    xwinR = XR(:,t0:n);
    xpsdL = abs(xwinL(2:end,:)).^2;
    xpsdR = abs(xwinR(2:end,:)).^2;
    xnt = size(xpsdL,2); %L/R same size
    for i = 1:xnt
        cent(i,1) = sum(xpsdL(:,i).*f(2:end))./sum(xpsdL(:,i));
        cent(i,2) = sum(xpsdR(:,i).*f(2:end))./sum(xpsdR(:,i));
    end
    E_fastL = sum(abs(xwinL).^2) + sum(abs(xwinL(2:end-1,:)).^2);
    E_fastR = sum(abs(xwinR).^2) + sum(abs(xwinR(2:end-1,:)).^2);
    E_slow(1) = max(E_fastL);
    E_slow(2) = max(E_fastR);
    E_normL = E_fastL./E_slow(1);
    E_normR = E_fastR./E_slow(2);
    A_normL = sqrt(E_normL);
    A_normR = sqrt(E_normR);
        
    if cent(end,1) < Fctf
        k(1) = exp(-(A_normL(end)-min(A_normL))/tau)*kv;
    else
        k(1) = kc;
    end
    
    if cent(end,2) < Fctf
        k(2) = exp(-(A_normR(end)-min(A_normR))/tau)*kv;
    else
        k(2) = kc;
    end
        
    Gexp(1) = A_normL(end).^(k(1));
    Gexp(2) = A_normR(end).^(k(2));
    %%%%
    
    %%run compression%%
    t0 = max(1,n-nbin_cmp+1);
    xwinL = XL(:,t0:n);
    xwinR = XR(:,t0:n);
    xpsdL = abs(xwinL(2:end,:)).^2;
    xpsdR = abs(xwinR(2:end,:)).^2;
    xnt = size(xpsdL,2); %L/R same size

    E_fastL = sum(abs(xwinL).^2) + sum(abs(xwinL(2:end-1,:)).^2);
    E_fastR = sum(abs(xwinR).^2) + sum(abs(xwinR(2:end-1,:)).^2);
    E_slow(1) = max(E_fastL);
    E_slow(2) = max(E_fastR);
    E_normL = E_slow(1)./max(1e-10,E_fastL);
    E_normR = E_slow(2)./max(1e-10,E_fastR);
    A_normL = sqrt(E_normL);
    A_normR = sqrt(E_normR);
        
    Gcmp(1) = min(gain_bounds_A(2),max(gain_bounds_A(1), A_normL(end)));
    Gcmp(2) = min(gain_bounds_A(2),max(gain_bounds_A(1), A_normR(end)));
    %%%%
    
    G(1) = Gexp(1)*Gcmp(1);
    G(2) = Gexp(2)*Gcmp(2);
    GmatL = repmat(G(1),length(f),1);
    GmatR = repmat(G(2),length(f),1);
    if ear == 0
        xw_enhL = xwinL(:,end).*GmatL;
        xw_enhR = xwinR(:,end).*GmatR;
    elseif ear == 1
        xw_enhL = xwinL(:,end).*GmatL;
        xw_enhR = xwinR(:,end).*GmatL;
    elseif ear == 2
        xw_enhL = xwinL(:,end).*GmatR;
        xw_enhR = xwinR(:,end).*GmatR;
    end
    xt_enhL = real(ifft([xw_enhL; conj(flipud(xw_enhL(2:end-1,:)))])); %samples-by-time matrix
    xt_enhL = xt_enhL.*hann(NFFT);
    xt_enhR = real(ifft([xw_enhR; conj(flipud(xw_enhR(2:end-1,:)))])); %samples-by-time matrix
    xt_enhR = xt_enhR.*hann(NFFT);
    xt0 = max(1,length(XT_enh_L)-NFFT/2+1); %L/R same length
    if n == 1
        XT_enh_L = cat(1,XT_enh_L,zeros(NFFT,1));
        XT_enh_R = cat(1,XT_enh_R,zeros(NFFT,1));
    else
        XT_enh_L = cat(1,XT_enh_L,zeros(NFFT/2,1));
        XT_enh_R = cat(1,XT_enh_R,zeros(NFFT/2,1));
    end
    XT_enh_L(xt0:end) = XT_enh_L(xt0:end)+xt_enhL;
    XT_enh_R(xt0:end) = XT_enh_R(xt0:end)+xt_enhR;
end

y(:,1) = XT_enh_L;
y(:,2) = XT_enh_R;
end