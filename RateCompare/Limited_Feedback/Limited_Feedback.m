function SE = Limited_Feedback(H_total, Chan, OFDM, BS, UE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this algorithm from the paper:
% Limited Feedback Hybrid Precoding for Multi-User Millimeter Wave Systems
% Ahmed Alkhateeb, Geert Leus, and RobertW. Heath, Jr.,
% the authors have assumed analog beamforming at the receive side. We keep
% the same assumption also here also in the other baselines we use multiple
% RF chains at the UE side
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
% For first way of implementation
% WRF = zeros(Nr,U);
% FRF = zeros(Nt,U);

% For second way of implementation 
WRF_ = zeros(Nr,U);
FRF_ = zeros(Nt,U);
FBB = zeros(U,U,K);

D_t = zeros(Nt, Nt);
delta_t = 2*pi/Nt;
for j=1:Nt
    D_t(:,j) = (1/sqrt(Nt))*exp(1i*(0:Nt-1)*(j-1)*delta_t);
end

D_r = zeros(Nr, Nr);
delta_r = 2*pi/Nr;
for j=1:Nr
    D_r(:,j) = (1/sqrt(Nr))*exp(1i*(0:Nr-1)*(j-1)*delta_r);
end

% average the channel over the subcarriers
H_mean = squeeze(mean(H_total,3));
% gain_max = 0;
% for u=1:U
%     H_u = H_mean(:,:,u);
%     for i=1:Nt
%         for j=1:Nr
%             gain = abs(D_r(:,j)'*H_u*D_t(:,i));
%             if gain>=gain_max
%                 gain_max = gain;
%                 WRF(:,u) = D_r(:,j);
%                 FRF(:,u) = D_t(:,i);
%             end
%         end
%     end
%     gain_max = 0;
% end

% Finding FRF, WRF
for u=1:U
    H_u = H_mean(:,:,u);
    gain = abs(D_r'*H_u*D_t);
    [idx, idy]=find(gain==max(gain(:)));
    WRF_(:,u) = D_r(:,idx);
    FRF_(:,u) = D_t(:,idy);
end

% Finding FBB
for k=1:K
    HH = [];
    for u=1:U
        hbar_u = (WRF_(:,u)'*H_total(:,:,k,u)*FRF_)';
        HH = [HH, hbar_u];
    end
    HH = HH';
    FBB(:,:,k) = HH'*pinv(HH*HH');
end

% Normalization
for k=1:K
    for u=1:U
        FBB(:,u,k) = (FBB(:,u,k))/(norm(FRF_*FBB(:,u,k), 'fro'));
    end
end

% Spectral efficiency calculation
SE = zeros(1, U);

% Noise power calculation [Fully connected Architecture]
Pn_Full = LimFeed_NoisePowerCal(Chan, H_total, OFDM, FRF_, FBB, WRF_);

for i = 1:U
    for k=1:K
        Pu_Full = WRF_(:,i)'*H_total(:,:,k,i)*FRF_*...
                  FBB(:,i,k)*FBB(:,i,k)'*...
                  FRF_'*H_total(:,:,k,i)'*WRF_(:,i);
    
        SE(1,i) = SE(1,i) + ...
            log2(1 + pinv(Pn_Full(k,i))*Pu_Full)/OFDM.nfft;
    end
end




end