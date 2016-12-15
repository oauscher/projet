function out = noise_awgn(Eb_N0_dB,rr,signalComplex,upSamplingFactor)
    SNR = (10^(Eb_N0_dB(rr)/10)); %EB/No
    nt = 1/sqrt(2)*[randn(1,length(signalComplex)) + 1i*randn(1,length(signalComplex))]; 
    Esym=sum(abs(signalComplex).^2)/(length(signalComplex));
    N_0 = Esym/SNR;
    noise = sqrt(N_0)* nt;
    out = 1/sqrt(upSamplingFactor)*signalComplex + noise; %%%%%adding awgn