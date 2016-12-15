%%%%%%%MSK like receiver, we take only c0 into account ( and not c1)%%%%%
function outPreco=viterbiC1Preco(r0n) %%%%%mainly a decoder%%%%%

outPreco=zeros(1,length(r0n));
outPreco(1:2:end)=-1i*r0n(1:2:end); %%%odd
outPreco(2:2:end)=-r0n(2:2:end); %%%even

outPreco = sign(real(outPreco));