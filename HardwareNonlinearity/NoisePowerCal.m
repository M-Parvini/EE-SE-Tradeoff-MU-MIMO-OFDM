function Pn = NoisePowerCal(Chan, OFDM, FRF, FBB, WRF, WBB, H, alpha, ...
    rho_d, rho_a, beta, Ruu, Rss, BO)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of RF chains
UENRF = OFDM.UERFchain;
BSNRF = OFDM.BSRFchain;
% Number of OFDM symbols
Ns = OFDM.numStreams;
% Number of UEs
U = OFDM.nUEs;
% noise
sigma = Chan.NoisePower;
% OFDM subcarrier numbers
K = OFDM.nfft;
% Pin
Pin = OFDM.Pin;
%  kappas
kappa1 = alpha*(1-rho_d)*(1-rho_a)/sqrt(BO);
kappa2 = alpha*(1-rho_a)/sqrt(BO);
kappa3 = (1-rho_a);

% total noise power
Pn = zeros(Ns,Ns,K,U);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inter UE interference
Pn1 = zeros(Ns,Ns,K,U);
for k=1:K
    for j=1:U
        for i=1:U
            if j~=i
                Pn1(:,:,k,j) = Pn1(:,:,k,j) + Pin/(U*Ns)*kappa1^2*...
                        WBB(:,:,k,j)'*WRF(:,:,j)'*H(:,:,k,j)*...
                        FRF*FBB(:,:,k,i)*FBB(:,:,k,i)'*...
                        FRF'*H(:,:,k,j)'*...
                        WRF(:,:,j)*WBB(:,:,k,j);
            end
        end
    end
end

% AWGN Noise
Pn2 = zeros(Ns,Ns,K,U);
for u=1:U
    for k=1:K
        Pn2(:,:,k,u) = sigma*kappa3^2*...
            WBB(:,:,k,u)'*(WRF(:,:,u)'*WRF(:,:,u))*WBB(:,:,k,u);
    end
end

% DAC Noise
Pn3 = zeros(Ns,Ns,K,U);
for u=1:U
    for k=1:K
        Pn3(:,:,k,u) = kappa2^2*rho_d*(1-rho_d)*WBB(:,:,k,u)'*...
        WRF(:,:,u)'*H(:,:,k,j)*FRF*Ruu(:,:,k)*FRF'*H(:,:,k,j)'*...
        WRF(:,:,u)*WBB(:,:,k,u);
    end
end

% PA noise
Pn4 = zeros(Ns,Ns,K,U);
for u=1:U
    for k=1:K
        Pn4(:,:,k,u) = kappa3^2*((beta-alpha^2)/BO)*WBB(:,:,k,u)'*...
            WRF(:,:,u)'*H(:,:,k,j)*Rss(:,:,k)*H(:,:,k,j)'*WRF(:,:,u)*WBB(:,:,k,u);
    end
end

% ADC noise
Pn5 = zeros(Ns,Ns,K,U);
for u=1:U
    for k=1:K
        Pn5(:,:,k) = (rho_a*(1-rho_a)*(beta)/BO)*...
            WBB(:,:,k,u)'*WRF(:,:,u)'*H(:,:,k,j)*Rss(:,:,k)*...
            H(:,:,k,j)'*WRF(:,:,u)*WBB(:,:,k,u)...
            +sigma*rho_a*(1-rho_a)*WBB(:,:,k,u)'*(WRF(:,:,u)'...
            *WRF(:,:,u))*WBB(:,:,k,u);
    end
end

Pn = Pn1+Pn2+Pn3+Pn4+Pn5;
end
