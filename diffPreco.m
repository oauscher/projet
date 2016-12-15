function nrzPreco=diffPreco(nrz) 
%     nrzPreco=zeros(1,length(nrz));
%     nrzPreco(3:2:end-1)= -mod((nrz(3:2:end-1)-nrz(2:2:end-2) + 1),4);
%     nrzPreco(2:2:end)= mod((nrz(2:2:end)-nrz(1:2:end-1) + 1),4);
%     nrzPreco(1)=nrz(1);
%     nrzPreco(find(nrzPreco==3))=-1;
%     nrzPreco(find(nrzPreco==-3))=1;
      nrzPreco=zeros(1,length(nrz));
      nrzPreco(1)=nrz(1);
      nrzPreco(2:1:end)=nrz(2:1:end).*nrz(1:1:end-1);
      for i = 2:length(nrzPreco)
        nrzPreco(i)=(-1)^i*nrzPreco(i);
      end