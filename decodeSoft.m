function [ dv, s ] = decodeSoft( r, H, stdNoise, it )
%   
%   GOES WITH MAIN3
%
%   Soft decoding of LDPC codeword
%
%   r  :  received codeword
%   H  :  parity-check matrix
%   it :  number of iteration
%
%   dv :  decoded vector
%   s  :  syndrom

%%%%%%%%%%%%%%%% Developer's Corner
%%%   1- pb dans Chouinard : son Qij0 est en fait Qij1 et vice-versa
%%%      vérifié avec original decodeprobdomain from FE qui corrige bien
%%%  2- Horizontal step correcte  21/11
%%%   3- Vertical step added 21/11
%%%   666- Everything seems ok 21/11  Further tests needed
%%%  
%%%%%%%%%%%%%%%%

%% Initialization
[k,n] = size(H);
a=2;
Pi1 =  1./(1+exp(-2*a*r./stdNoise));
Pi0 = 1-Pi1;

Pij1 = H.*repmat(Pi1,k,1);
Pij0 = H.*repmat(Pi0,k,1);
Qij0 = zeros(k,n);
Qij1 = zeros(k,n);

% ------------ Loop ------------
for iter = 1:it
    %% Horizontal Step
    % SLIDE 33 CHOUINARD
    % For each row
    for i=1:k
        % Find non-zeros elements in H not to get a product of 0 in deltaQij
        indC = find(H(i,:)~=0);
        % For each non-zero element
        for l = 1:length(indC)
            deltaQij = 1; % Initialize deltaQij coeff for each iter
            % Get the product on same row of the deltaPij matrix to fill deltaQij
            for m = 1:length(indC)
                % Prevent product by current coeff
                if m~=l
                    deltaPij = Pij0(i,indC(m)) - Pij1(i,indC(m));
                    deltaQij = deltaQij * deltaPij;
                end
            end
            Qij0(i,indC(l)) = 0.5*(1+deltaQij);
            Qij1(i,indC(l)) = 0.5*(1-deltaQij);
        end
    end % for each row
      
    %% Vertical Step
    % SLIDE 34 CHOUINARD
    % For each column
    for j=1:n
        % Find non-zeros elements in H not to get a product of 0 in
        indL = find(H(:,j)~=0);
        % For each non-zero element in the j column
        for l = 1:length(indL)
            % Get the product on same column of the Qij matrix
            for m = 1:length(indL)
                xQij0 = 1;
                xQij1 = 1;
                % Prevent product by current coeff
                if m~=l
                    xQij0 = xQij0 * Qij0(indL(m),j);
                    xQij1 = xQij1 * Qij1(indL(m),j);
                end
            end
            % Temporary coefficients to compute normalized probability matrix
            tmp0 = Pi0(j) * xQij0;
            tmp1 = Pi1(j) * xQij1;
            
            % New probability matrix
            Pij0(indL(l),j) = tmp0 / (tmp0 + tmp1);
            Pij1(indL(l),j) = tmp1 / (tmp0 + tmp1);
        end
        % Temporary coefficients to compute normalized probability vectors
        tmpPi0 = Pi0(j)*prod(Qij0(indL,j));
        tmpPi1 = Pi1(j)*prod(Qij1(indL,j));
        
        %New probability vectors
        Pi0(j) = tmpPi0 / (tmpPi0 + tmpPi1);
        Pi1(j) = tmpPi1 / (tmpPi0 + tmpPi1);
        
        % LLR for decision
        llr = log(Pi1./Pi0);
        cEsti = (sign(llr)+1)/2;
    end % for each column 
end % iterations
% ------------ Loop End ------------%

%% Finalization
%Decoded vector & Syndrom
dv = cEsti(1:k);
s = mod(cEsti*H',2); 
% llr
% Pi1
% Pi0
end % function


