function r0n=pulsShapeReceiver(c,out,upSamplingFactor) %%%%%%Conv with c0
    y=conv(c,out);
    r0n=y(length(c):upSamplingFactor:end); %Nco = length co, removing it and sampling at frequency Nb (= Ts)