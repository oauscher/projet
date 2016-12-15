clc;
clear all;
%%%%%%%%%%%%%%%%%modulation parameters%%%%%%%%%%%%

h=0.5; %modulation indice
Ts=1;
upSamplingFactor=32;
Tp=Ts/upSamplingFactor;% sampling period
BT=0.3;
L=4;
sig=sqrt(log(2))/(2*pi*BT);

t=0:Tp:L*Ts-Tp;%time basis for to generate the waveform
tpp=0:Tp:2*L*Ts-Tp; %time basis for the c0 approx
%%%%%%%%%%%%impulse response of the gaussian filter%%%%%%%%%%%%%%%%

y=0.5*(erf((1/(sqrt(2)*sig))*((t-3*Ts/2)/Ts))-erf((1/(sqrt(2)*sig))*((t-5*Ts/2)/Ts)));
x=rectpuls(t-2*Ts,4*Ts);
g=y.*x;
g=g;
q = cumsum(g)/sum(g);%%% normalisation to have int q = 1

phi=q; 
for ii=1:length(q)
     phi(ii+length(q))= phi(ii+length(q)-1)-(g(ii)/sum(g));
end
% 
p0 = sin(pi*h*phi)/sin(pi*h);
% %p1 = p(t+Ts)
p1 = [p0 zeros(1, upSamplingFactor)];
p1 = p1(upSamplingFactor+1:end);
% 
p2 = [p0 zeros(1, 2*upSamplingFactor)];
p2 = p2(2*upSamplingFactor+1:end);
% 
p3 = [p0 zeros(1, 3*upSamplingFactor)];
p3 = p3(3*upSamplingFactor+1:end);

p4 = [p0 zeros(1, 4*upSamplingFactor)];
p4 = p4(4*upSamplingFactor+1:end);

p5 = [p0 zeros(1, 5*upSamplingFactor)];
p5 = p5(5*upSamplingFactor+1:end);
% 
% 
ct=p0.*p1.*p2.*p3;
c=fliplr(ct); %c(-t)

ct1=p0.*p2.*p5.*p3;
c1=fliplr(ct1); %c(-t)


%%%%%%%inputBit%%%%%%%%%%%%%%%%%%%%%%%


Lb = 9600+4; % Lb bits %%% multiple of 12 (+4 because for bits are deleted) size of the input bit sequence
Eb_N0_dB = 1:10;

for rr = 1:length(Eb_N0_dB)
H = comm.PNSequence('SamplesPerFrame',Lb);
inputBits = H.step';

inputBitsForTheEnd = inputBits; %%%%%%used for performance computation

encodedBits = encode( inputBits(4:end-1), 6, 12 );
encodedBits2 = [inputBits(1:3) encodedBits inputBits(length(inputBits))];
inputBits = encodedBits2; %%%% pour que le decodeur viterbi est la bonne taille de sequence

%inputBits =  [1 0 1 0 1 1 0 1 0 1 1 1 1 0 1 0 1 1 1 1 ] ;
nrz = 2*encodedBits2 - 1; %%%%%%%


nrzPreco = diffPreco(nrz);

signalComplex = pulsShap(nrz,upSamplingFactor,g,h);
signalComplexPreCo = pulsShap(nrzPreco,upSamplingFactor,g,h);

%%%%%%%%spectrumù%%%%%%%%%%%%%%%%%%%%%

%  L=length(signalSent);
%  NFFT = 1024;
%  X = abs(fftshift(fft(signalSent,NFFT)));
%  f=-1/2:1/(NFFT-1):1/2;

%%%%%%%%%%%%%%%%%%%%%canal%%%%%%%%%%%%%%%%

out = noise_awgn(Eb_N0_dB,rr,signalComplex,upSamplingFactor);
outPreCo = noise_awgn(Eb_N0_dB,rr,signalComplexPreCo,upSamplingFactor);
varPreCo = puissanceBruitBlanc(Eb_N0_dB,rr,signalComplexPreCo,upSamplingFactor);
%%%%%%%%%%%%%%%%%%%%%%%%%receiver%%%%%%%%%%%%%%%%

r0n=pulsShapeReceiver(c,out,upSamplingFactor);
r1n=pulsShapeReceiver(c1,out,upSamplingFactor);

r0nPreCo=pulsShapeReceiver(c,outPreCo,upSamplingFactor);
r1nPreCo=pulsShapeReceiver(c1,outPreCo,upSamplingFactor);

%%%%%%%%%%%%%%%viterbi decoder%%%%%%%%%%%%%%%%%%%%%

out=viterbi(inputBits,r0n,r1n);
out=diffDeco(out);

out2=viterbiC1(r0n);
outPreCoC1=viterbiWhenPreco(inputBits,r0nPreCo,r1nPreCo);
outPreCo=viterbiC1Preco(r0nPreCo);
% nrz(4:end-1)
nrz=nrz(4:end-1);

out=out(3:end-1);
out2=out2(3:end-1);
outPreCoC1=outPreCoC1(4:end-1);
outPreCo=outPreCo(4:end-1);

%%%%decodage de canal%%%%%
H = [ 1 1 0 1 0 0 0 0 1 1 0 1;
      1 1 1 1 1 0 0 1 0 0 0 0;
      1 0 1 0 1 1 0 0 0 1 1 0;
      0 0 0 1 1 1 1 0 1 0 1 0;
      0 1 1 0 0 0 1 1 0 0 1 1;
      0 0 0 0 0 1 1 1 1 1 0 1];

NrzEncodedBits = inputBitsForTheEnd(4:end-1);  
decodedSeq = decodeLDPC( outPreCoC1 ,H,varPreCo, 10 ); %%%last arg, iteration number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CODAGE DE CANAL AVEC TEB5%%%%%%%%%%%%%%%%%%%%%


 if rr == 1 %initialisation du vecteur
        TEB = (length(out)-length(find(out == nrz)))/length(nrz);
        TEB2 = (length(out2)-length(find(out2 == nrz)))/length(nrz);
        TEB3 = (length(outPreCo)-length(find(outPreCo == nrz)))/length(nrz);
        TEB4 = (length(outPreCoC1)-length(find(outPreCoC1 == nrz)))/length(nrz);
        TEB5 =  (length(decodedSeq)-length(find(decodedSeq == NrzEncodedBits)))/length(NrzEncodedBits);
    else
        TEB = [TEB, (length(out)-length(find(out == nrz)))/length(nrz)];
        TEB2 = [TEB2, (length(out2)-length(find(out2 == nrz)))/length(nrz)];
        TEB3 = [TEB3, (length(outPreCo)-length(find(outPreCo == nrz)))/length(nrz)];
        TEB4 = [TEB4, (length(outPreCoC1)-length(find(outPreCoC1 == nrz)))/length(nrz)];
        TEB5 = [TEB5, (length(decodedSeq)-length(find(decodedSeq == NrzEncodedBits)))/length(NrzEncodedBits)];
 end
end
 
EbNoLin=10.^(Eb_N0_dB/10);
gg=erfc(sqrt(EbNoLin));%theo curve
gg2=0.5*erfc(sqrt(EbNoLin));%theo curve

figure(1)
semilogy(Eb_N0_dB,gg);
hold on
semilogy(Eb_N0_dB,gg2);
hold on
semilogy(Eb_N0_dB,TEB,'r');
hold on 
semilogy(Eb_N0_dB,TEB2,'b');
hold on 
semilogy(Eb_N0_dB,TEB3,'g');
hold on
semilogy(Eb_N0_dB,TEB4);
hold on 
semilogy(Eb_N0_dB,TEB5);
grid on
legend('erfc(sqrt(Eb/No))', '0.5.erfc(sqrt(Eb/No))','2 PAM approx','1 PAM approx','with Precoding 1 PAM','with Precoding 2 PAM','with Precoding 2 PAM and LDPC');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for GMSK modulation');
