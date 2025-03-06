function results = MIMO_OFDM(OFDM, Chan, BS, UE, SNR)
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
OFDM.SNR = db2pow(SNR);
Chan.NoisePower = 1/OFDM.SNR;
H_total = zeros(Nr,Nt,K,U);

%% MIMO channel modeling for UEs
for u=1:U
    [Chan, H_fc,Atfc,Arfc,gains] = MassiveMimoChannel(Chan, OFDM, BS, UE);
    H_total(:,:,:,u) = H_fc;
    Atfc_total(:,:,u) = Atfc;
    Arfc_total(:,:,u) = Arfc;
    gains_total(:,:,u) = gains;
end

%% Fully Digital
% [Dig_SE_indep] = Fully_Digital_Independent(H_total, Chan, OFDM, BS, UE);
% [Dig_SE_coop] = Fully_Digital_Cooperative(H_total, Chan, OFDM, BS, UE);
[DigRCD_SE] = Fully_Digital_RCD(H_total, Chan, OFDM, BS, UE);

%% Proposed solution [RCD-based]
[Pro_SE] = Hybrid_Proposed(H_total, Chan, OFDM, BS, UE);

%% Khalid --> approximated RCD
[RCD_SE] = KhalidRCD(H_total, Chan, OFDM, BS, UE,Atfc_total,Arfc_total,...
    gains_total);

%% MO Hybrid Precoder Design
[MO_SE] = MO_Hybrid(H_total, Chan, OFDM, BS, UE);

%% BD Hybrid Precoder Design
[BD_SE] = BD_Hybrid(H_total, Chan, OFDM, BS, UE);

%% Limited_Feedback
[LimFeed_SE] = Limited_Feedback(H_total, Chan, OFDM, BS, UE);

%%% Saving the results
results.Dig_SE_RCD = DigRCD_SE;
results.Proposed_SE = Pro_SE;
results.KhalidRCD_SE = RCD_SE;
results.Lim_Feed = LimFeed_SE;
results.BD_SE  = BD_SE;
results.MO_SE  = MO_SE;
end