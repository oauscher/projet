%%%%%%%MSK like receiver, we take only c0 into account ( and not c1)%%%%%
function out2=viterbiC1(r0n) %%%%%mainly a decoder%%%%%

a2n1=r0n(2:1:end);

a2n=r0n(1:1:end-1);

out2=-1i*a2n1.*a2n;

out2(1:2:end)=-out2(1:2:end);

out2 = sign(real(out2));