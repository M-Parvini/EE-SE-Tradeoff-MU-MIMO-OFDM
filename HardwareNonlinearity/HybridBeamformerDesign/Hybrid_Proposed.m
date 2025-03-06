function [FRF,FBB,WRF,WBB] = Hybrid_Proposed(H_total, Chan, OFDM, BS, UE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% Number of RF chains
Nr_rf = OFDM.UERFchain;
Nt_rf = OFDM.BSRFchain;
% Number of OFDM symbols
Ns = OFDM.numStreams;
% Number of UEs
U = OFDM.nUEs;
% noise
sigma = Chan.NoisePower;
% central frequency fc
fc = Chan.fc;
% bandwidth
fs = OFDM.BW;
% max time delay
DelaySpread = Chan.delay_spread;
% OFDM subcarrier numbers
K = OFDM.nfft;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of precoder and combiner matrices
FRF = zeros(Nt,Nt_rf);
FBB = zeros(Nt_rf,Ns,K,U);
WRF = zeros(Nr,Nr_rf,U);
WBB = zeros(Nr_rf,Ns,K,U);
H_G = zeros(Nt,U*Nr_rf,K);
% Start of algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RF Combiner Design
D_r = zeros(Nr, Nr);
delta_r = 2*pi/Nr;
for j=1:Nr
    D_r(:,j) = (1/sqrt(Nr))*exp(1i*(0:Nr-1)*(j-1)*delta_r);
end

for u=1:U
    temp = 0;
    for k=1:K
        Hk = H_total(:,:,k,u);
        temp = temp + sum(abs(D_r'*Hk),2);
    end
    [~, idx] = sort(temp, 'descend');
    WRF(:,:,u) = D_r(:,idx(1:Nr_rf));

    for k=1:K
        H_u = H_total(:,:,k,u);
        H_G(:,(u-1)*Ns+1:u*Ns,k) = (WRF(:,:,u)'*H_u).';
    end
end

gain_old = 0;
for k=1:K
    gain_new = 0;
    H_Gk = H_G(:,:,k).';
    coeff = 1/sqrt(Nt);
    FRF_ = coeff*exp(1i*angle(H_Gk'));
    for u=1:U
        gain_new = gain_new + norm(WRF(:,:,u)'*H_total(:,:,k,u)*FRF_,'fro');
    end
    if gain_new >= gain_old
        gain_old = gain_new;
        FRF = FRF_;
    end
end



% Finding the equivalent channel
Heq = zeros(Nr_rf,Nt_rf,K,U);
for u=1:U
    for k=1:K
        Heq(:,:,k,u) = WRF(:,:,u)'*H_total(:,:,k,u)*FRF;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%% Block DIagonalization %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Precoders and Combiners
for m = 1:K
    for j = 1:U
        H_tildeJ = [];
        for k = setdiff(1:U,j)
            H_tildeJ = [H_tildeJ(:,:);   Heq(:,:,m,k)]; 
        end

        Fa = pinv(H_tildeJ'*H_tildeJ + ((Nr_rf*sigma))*eye(Nt_rf));
        [U_,~,V_] = svd(Heq(:,:,m,j)*Fa);
        
        Fb = V_(:,1:Ns);
        FBB(:,:,m,j) = Fa*Fb;
        FBB(:,:,m,j) = FBB(:,:,m,j)/norm(FBB(:,:,m,j),'fro');
        WBB(:,:,m,j) = U_(:,1:Ns);
        WBB(:,:,m,j) = WBB(:,:,m,j)/norm(WBB(:,:,m,j),'fro');
    end
end

% Normalization
for k=1:K
    for u=1:U
        FBB(:,:,k,u) = FBB(:,:,k,u)/(norm(FRF*FBB(:,:,k,u),'fro'));
        WBB(:,:,k,u) = WBB(:,:,k,u)/(norm(WRF(:,:,u)*WBB(:,:,k,u),'fro'));
    end
end

end