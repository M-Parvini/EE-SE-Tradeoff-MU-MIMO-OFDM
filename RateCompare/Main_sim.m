clc
clear
% close all
rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);
tic
add_paths
%%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
NumSim = 200;
SNR = -20:5:10;                  % Different input power [dBm]
TransmitAntennas = 64;
ReceiveAntennas = 16;
numOFDMSym = 4;                % Number of OFDM symbols per UE (Ns)
UERFChain = 4;
numUE = 8;                     % Number of UEs
%%%%%%%%%%%%%%%% Simulation chain %%%%%%%%%%%%%%%%%%%%
% Fully connected
Results_DigRCD = zeros([NumSim, length(SNR), length(numUE)]);
Results_Proposed = zeros([NumSim, length(SNR), length(numUE)]);
Results_KhalidRCD = zeros([NumSim, length(SNR), length(numUE)]);
Results_BDSE = zeros([NumSim, length(SNR), length(numUE)]);
Results_MOSE = zeros([NumSim, length(SNR), length(numUE)]);
Results_LimFeedSE = zeros([NumSim, length(SNR), length(numUE)]);

% Loop
for nt = 1:length(TransmitAntennas)
    Nt = TransmitAntennas(nt);
    Nr = ReceiveAntennas;
    for u=1:length(numUE)
        BSRFChain = numUE(u)*UERFChain;   % RF chains (NRF) U*UERFChain
        for snr = 1:length(SNR)
            fprintf('::::::::::::::::::\n')
            fprintf(['SNR =', num2str(SNR(snr)), '\n'])
            [OFDMParams, ChanParams, BSParams, UEParams] = ...
                    InitializeParams(SNR, Nt, Nr, BSRFChain, UERFChain, ...
                    numOFDMSym, numUE(u));
            parfor SimId = 1:NumSim
                %%% Simulation Cycle
                results = ...
                    MIMO_OFDM(OFDMParams, ChanParams, BSParams, UEParams, SNR(snr));
    
                Results_DigRCD(SimId, snr, u) = sum(real(results.Dig_SE_RCD));
                Results_Proposed(SimId, snr, u) = sum(real(results.Proposed_SE));
                Results_KhalidRCD(SimId, snr, u) = sum(real(results.KhalidRCD_SE));
                Results_BDSE(SimId, snr, u) = sum(real(results.BD_SE));
                Results_MOSE(SimId, snr, u) = sum(real(results.MO_SE));
                Results_LimFeedSE(SimId, snr, u) = sum(real(results.Lim_Feed));
            end
        end
    end
end
save('TotalResults')


toc
% delete(poolobj)