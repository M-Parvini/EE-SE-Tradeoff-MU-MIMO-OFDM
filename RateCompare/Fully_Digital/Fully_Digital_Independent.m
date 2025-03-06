function SE = Fully_Digital_Independent(H_total, Chan, OFDM, BS, UE)
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

%%% Precoders and Combiners
WD = zeros(Nr,Ns,K,U);
FD = zeros(Nt,Ns,K,U);

for u=1:U
    for k=1:K
        [U_,~,V_] = svd(H_total(:,:,k,u));
        WD(:,:,k,u) = U_(:,1:Ns);
        FD(:,:,k,u) = V_(:,1:Ns);
    end
end

%% Capacity calculation
SE = zeros(1, U);

% Noise power calculation [Fully connected Architecture]
Pn_Full = Fully_Digital_NoisePowerCal(Chan, H_total, OFDM, FD, WD);

for i = 1:U
    for k=1:K
        Pu_Full = WD(:,:,k,i)'*H_total(:,:,k,i)*...
                  FD(:,:,k,i)*FD(:,:,k,i)'*...
                  H_total(:,:,k,i)'*WD(:,:,k,i);
    
        SE(1,i) = SE(1,i) + ...
            log2(det(eye(Ns) + pinv(Pn_Full(:,:,k,i))*Pu_Full))/OFDM.nfft;
    end
end

end