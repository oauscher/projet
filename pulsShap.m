function signalComplex = pulsShap(nrz,upSamplingFactor,g,h)
    nrz_dirac = upsample(nrz,upSamplingFactor);
    nrz_gauss = conv(g,nrz_dirac);
    nrz_gauss = nrz_gauss(1:1:length(nrz_gauss)-length(g)+1);%%% cut the end duz to the conv ope
    Pt = pi*h*cumsum(nrz_gauss)/sum(g);% integration & implementation of pi h
    signalComplex = exp(1i*Pt);