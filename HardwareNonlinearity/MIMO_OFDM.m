function results = MIMO_OFDM(OFDM, Chan, BS, UE, Pin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASSIVE MIMO
% Beam Squint Effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
Nt = BS.nAntenna;       % # transmit antennas
Nr = UE.nAntenna;       % # receive antennas
Ns = OFDM.numStreams;   % # transmit symbols
K  = OFDM.nfft;
U = OFDM.nUEs;
OFDM.Pin = db2pow(Pin-30);
BO = db2pow(OFDM.PABO);
H_total = zeros(Nr,Nt,K,U);

%% MIMO channel modeling for UEs
for u=1:U
    [Chan, H_fc] = MassiveMimoChannel(Chan, OFDM, BS, UE);
    H_total(:,:,:,u) = H_fc;
end

%% Proposed solution [RCD-based]
[FRF,FBB,WRF,WBB] = Hybrid_Proposed(H_total, Chan, OFDM, BS, UE);

% DAC
[rho_d, Ruu] = DACModel(OFDM, FBB);

% Power Amplifier
[Rss, beta, alpha] = PowerAmplifier(OFDM, FRF, Ruu, rho_d, Nt);

% ADC
[rho_a, ~] = ADCModel(OFDM);

%% Capacity calculation
results.SE = zeros(1, U);

% Noise power calculation [Fully connected Architecture]
Pn_Full = NoisePowerCal(Chan, OFDM, FRF, FBB, WRF, WBB, H_total, ...
                        alpha, rho_d, rho_a, beta, Ruu, Rss, BO);

for i = 1:U
    for k=1:K
        Pu_Full = (1/BO)*alpha^2*(1-rho_d)^2*(1-rho_a)^2*...
                        OFDM.Pin/(U*Ns)*...
                        WBB(:,:,k,i)'*WRF(:,:,i)'*H_total(:,:,k,i)*...
                        FRF*FBB(:,:,k,i)*FBB(:,:,k,i)'*...
                        FRF'*H_total(:,:,k,i)'*...
                        WRF(:,:,i)*WBB(:,:,k,i);

        results.SE(1,i) = results.SE(1,i) + ...
            log2(det(eye(Ns) + pinv(Pn_Full(:,:,k,i))*Pu_Full))/OFDM.nfft;
    end
end

% power consumption
[P_FC, ~] = PowerConsumptionModel(Chan, OFDM, BS, UE, OFDM.Pin);
results.EE = results.SE/P_FC;


end