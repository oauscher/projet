function out=viterbiWhenPreco(inputBits,r0n,r1n)
length1=0;

length2=0;

length3=0;

length4=0;

l1=[1 1]; %viterbi decoder state

l2=[-1 1];

l3=[1 -1];

l4=[-1 -1];

for n=1:length(inputBits)/2-1 %NN length of inputbitstream 

g1=l1;

g2=l2;

g3=l3;

g4=l4;

length11=-imag(r0n(2*n+1))+real(r1n(2*n+1))+length1;

length12=-imag(r0n(2*n+1))-real(r1n(2*n+1))+length2;

length21=-imag(r0n(2*n+1))-real(r1n(2*n+1))+length3;

length22=-imag(r0n(2*n+1))+real(r1n(2*n+1))+length4;

length31=+imag(r0n(2*n+1))-real(r1n(2*n+1))+length1;

length32=+imag(r0n(2*n+1))+real(r1n(2*n+1))+length2;

length41=+imag(r0n(2*n+1))+real(r1n(2*n+1))+length3;

length42=+imag(r0n(2*n+1))-real(r1n(2*n+1))+length4;

if length11>length12

length1=length11;

l1=[g1,1];

else length1=length12;

l1=[g2,1];

end

if length21>length22

length2=length21;

l2=[g3,1];

else length2=length22;

l2=[g4,1];

end

if length31>length32

length3=length31;

l3=[g1,-1];

else length3=length32;

l3=[g2,-1];

end

if length41>length42

length4=length41;

l4=[g3,-1];

else length4=length42;

l4=[g4,-1];

end

g1=l1;

g2=l2;

g3=l3;

g4=l4;

length11=-real(r0n(2*n+2))+imag(r1n(2*n+2))+length1;

length12=-real(r0n(2*n+2))-imag(r1n(2*n+2))+length2;

length21=-real(r0n(2*n+2))-imag(r1n(2*n+2))+length3;

length22=-real(r0n(2*n+2))+imag(r1n(2*n+2))+length4;

length31=+real(r0n(2*n+2))-imag(r1n(2*n+2))+length1;

length32=+real(r0n(2*n+2))+imag(r1n(2*n+2))+length2;

length41=+real(r0n(2*n+2))+imag(r1n(2*n+2))+length3;

length42=+real(r0n(2*n+2))-imag(r1n(2*n+2))+length4;

if length11>length12

length1=length11;

l1=[g1,1];

else length1=length12;

l1=[g2,1];

end

if length21>length22

length2=length21;

l2=[g3,1];

else length2=length22;

l2=[g4,1];

end

if length31>length32

length3=length31;

l3=[g1,-1];

else length3=length32;

l3=[g2,-1];

end

if length41>length42

length4=length41;

l4=[g3,-1];

else length4=length42;

l4=[g4,-1];

end



end

path_length=max(max(length1,length2),max(length3,length4));

if length1==path_length

out=l1;

elseif length2==path_length

out=l2;

elseif length3==path_length

out=l3;

elseif length4==path_length

out=l4;

end
out(1:2:end)=-out(1:2:end);