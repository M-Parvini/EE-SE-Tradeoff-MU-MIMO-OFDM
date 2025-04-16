function [P_tot, P_u] = PowerConsumptionModel(Chan, OFDM, BS, UE, Ps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% Number of RF chains
UENRF = OFDM.UERFchain;
BSNRF = OFDM.BSRFchain;
% Number of OFDM symbols
Ns = OFDM.numStreams;
% central frequency fc
fc = Chan.fc;
% bandwidth
fs = OFDM.BW;
% fs = 1e9;
% number of users
U = OFDM.nUEs;
% max time delay
DelaySpread = Chan.delay_spread;
% OFDM subcarrier numbers
N = OFDM.nfft;
% DAC resolution
b = OFDM.DAClevels;
% PA PMAX
Pmax = db2pow(OFDM.Pmax-30);
% PA Gain
PAgain = db2pow(OFDM.PAGain);
% Input back-off
BO = OFDM.PABO;
epsilon = 1/(10^(BO/10));
% PA Doherty level
l = OFDM.Dohertylevel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% elements power
PBB = 200e-3;                   % Baseband unit power
PM = 0.3e-3;                     % Mixer
PLO = 22.5e-3;                  % Local Oscillator
PLP = 14e-3;                    % Low-path filter
PPS = 20e-3;                    % Phase shifter power
PDAC = 1.5e-5*2^b+9e-12*b*fs;   % DAC power consumption
% PADC = 1.8853e-15*2^b*fs*sqrt(1+(fs/560e6)^2);
PADC = 494e-15*2^b*fs;
PRFTx = 2*(PLP+PM+PDAC);
PRFRx = 2*(PLP+PM+PADC);
PPA = Power_Amplifier_power1(Pmax, epsilon, Nt, l);

%%% calculating the consumed powers (Tx+Rx-->*2)
P_tot   = PBB + PLO + BSNRF*PRFTx + Nt*PPA + BSNRF*Nt*PPS + Ps + ...
          PBB + PLO + U*(UENRF*PRFRx +          UENRF*Nr*PPS);

P_u   =  (PBB + PLO + BSNRF*PRFTx + Nt*PPA + BSNRF*Nt*PPS + Ps)/U + ...
          PBB + PLO + UENRF*PRFRx +          UENRF*Nr*PPS;

end

function PPA = Power_Amplifier_power1(Pmax, epsilon, Nt, l)
    
    if epsilon <= 1/(l^2)
        PPA = (4/(l*pi))*Pmax*sqrt(epsilon);
    else
        PPA = (4/(l*pi))*Pmax*((l+1)*sqrt(epsilon)-1);
    end

end

function PPA = Power_Amplifier_power2(Pmax, PAgain, Nt)
    Pout = PAgain;
    if Pout/Nt <= 0.25*Pmax
        PPA = (2/pi)*sqrt((Pout/Nt)*Pmax);
    else
        PPA = (6/pi)*sqrt((Pout/Nt)*Pmax) - 2*Pmax/pi;
    end

end
