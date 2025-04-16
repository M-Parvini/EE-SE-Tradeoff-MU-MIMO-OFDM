clc
clear
% close all
rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);
tic
add_paths
%%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
NumSim = 200;
Pin = -40:10:20;                   % Different input power [dBm]
TransmitAntennas = [64, 128];
ReceiveAntennas = 4;
numOFDMSym = 4;                    % Number of OFDM symbols per UE (Ns)
UERFChain = 4;
numUE = 2:10;                      % Number of UEs
% Hardware parameters
DAClevel = 2:5;                    % resolution (bits)
PABO = 0:0.5:8;                    % Power amplifier back-off
%%%%%%%%%%%%%%%% Simulation chain %%%%%%%%%%%%%%%%%%%%
ResultsSE_sum = zeros([NumSim, length(Pin), length(numUE), length(PABO), ...
    length(DAClevel), length(TransmitAntennas)]);

ResultsEE_sum = zeros([NumSim, length(Pin), length(numUE), length(PABO), ...
    length(DAClevel), length(TransmitAntennas)]);

ResultsSE_mean = zeros([NumSim, length(Pin), length(numUE), length(PABO), ...
    length(DAClevel), length(TransmitAntennas)]);

ResultsEE_mean = zeros([NumSim, length(Pin), length(numUE), length(PABO), ...
    length(DAClevel), length(TransmitAntennas)]);

% Loop
for nt = 1:length(TransmitAntennas)
    Nt = TransmitAntennas(nt);
    Nr = ReceiveAntennas;
    for u=1:length(numUE)
        BSRFChain = numUE(u)*UERFChain;   % RF chains (NRF) U*UERFChain
        for IBO = 1:length(PABO)
            for b = 1:length(DAClevel)
                for p = 1:length(Pin)
                    fprintf('::::::::::::::::::\n')
                    fprintf(['Pin =', num2str(Pin(p)), '\n'])
                    [OFDMParams, ChanParams, BSParams, UEParams] = ...
                            InitializeParams(Pin, Nt, Nr, BSRFChain, ...
                            UERFChain, numOFDMSym, numUE(u), ...
                            DAClevel(b), PABO(IBO));
                    for SimId = 1:NumSim
                        %%% Simulation Cycle
                        results = ...
                            MIMO_OFDM(OFDMParams, ChanParams, BSParams,...
                            UEParams, Pin(p));
            
                        ResultsSE_sum(SimId,p,u,IBO,b,nt) = sum(real(results.SE));
                        ResultsSE_mean(SimId,p,u,IBO,b,nt) = mean(real(results.SE));
                        ResultsEE_sum(SimId,p,u,IBO,b,nt) = real(results.EEtot);
                        ResultsEE_mean(SimId,p,u,IBO,b,nt) = mean(real(results.EEu));
                    end
                end
            end
        end
    end

end
save('TotalResults')


toc
% delete(poolobj)