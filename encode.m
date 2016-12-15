function [ encodedbits ] = encode( bits, k, n )
%k=6    n=12
G = [1 0 0 0 0 0 0 1 1 0 1 0;
    0 1 0 0 0 0 0 1 1 1 1 1;
    0 0 1 0 0 0 1 1 0 0 1 0;
    0 0 0 1 0 0 0 1 0 1 1 0;
    0 0 0 0 1 0 1 1 1 0 1 1;
    0 0 0 0 0 1 1 0 0 1 0 1];

bitscol = reshape(bits, k, length(bits)/k);
encodedbits = zeros(1,(length(bits)/k)*n);
for i = 1:size(bitscol,2)
    encodedbits(1,12*(i-1)+1:12*i) = mod(bitscol(:,i)'*G,2);
end
% matrice
% encodedbits = reshape(encodedbits,12,(length(bits)/k)*n);
end