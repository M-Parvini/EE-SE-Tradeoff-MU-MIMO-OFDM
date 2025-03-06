function [SE] = BD_Hybrid(H_total, Chan, OFDM, BS, UE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% Number of RF chains
NRF = OFDM.UERFchain;
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

[BD_WRF] = BD_Combiner_RF(H_total, Chan, OFDM, BS, UE);
[BD_FRF] = BD_Precoder_RF(H_total, Chan, OFDM, BS, UE, BD_WRF);
[BD_FBB, BD_WBB] = BD_Baseband(H_total, Chan, OFDM, BS, UE, BD_FRF, BD_WRF);

%%%
FRF = BD_FRF;
WRF = BD_WRF;
FBB = BD_FBB;
WBB = BD_WBB;

%% Capacity calculation
SE = zeros(1, U);

% Noise power calculation [Fully connected Architecture]
Pn_Full = BD_NoisePowerCal(Chan, H_total, OFDM, FRF, FBB, WRF, WBB);

for i = 1:U
    for k=1:K
        % Pu_Full = (OFDM.Pin/(U*Ns))*...
        %                 WBB(:,:,k,i)'*WRF(:,:,i)'*H_total(:,:,k,i)*FRF*...
        %                 FBB(:,:,k,i)*FBB(:,:,k,i)'*...
        %                 FRF'*H_total(:,:,k,i)'*WRF(:,:,i)*WBB(:,:,k,i);

        Pu_Full = WBB(:,:,k,i)'*WRF(:,:,i)'*H_total(:,:,k,i)*FRF*...
                  FBB(:,:,k,i)*FBB(:,:,k,i)'*...
                  FRF'*H_total(:,:,k,i)'*WRF(:,:,i)*WBB(:,:,k,i);
    
        SE(1,i) = SE(1,i) + ...
            log2(det(eye(Ns) + pinv(Pn_Full(:,:,k,i))*Pu_Full))/OFDM.nfft;
    end
end
end