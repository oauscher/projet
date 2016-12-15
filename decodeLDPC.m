function [ decodedBits ] = decodeLDPC( encodedBits ,H,stdNoise, it )

k=6;n=12;
encodedBitsCols = reshape(encodedBits,12,length(encodedBits)/12);
decodedBits = zeros(1,(length(encodedBits)*(k/n)));

for i=1:length(encodedBits)/12
    [dv,s] = decodeSoft(encodedBitsCols(:,i)',H,stdNoise, it );
    dv=dv(1:k);
    decodedBits(1,6*(i-1)+1:6*i) = dv;
end