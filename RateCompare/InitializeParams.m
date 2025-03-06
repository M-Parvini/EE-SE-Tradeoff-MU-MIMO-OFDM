function [OFDM, Chan, BS, UE] = InitializeParams(SNR, Nt, Nr, ...
    BSRFChain, UERFChain, numOFDMSym, numUE)

%%% Modulation parameters
OFDM.SNR = 0;                   % Input power
OFDM.SNRList = SNR;             % Input power list
OFDM.numStreams = numOFDMSym;   % Number of parallel data streams (Ns)
OFDM.BSRFchain = BSRFChain;     % Number of RF Chains (NRF)
OFDM.UERFchain = UERFChain;     % Number of RF Chains (NRF)
OFDM.nUEs = numUE;              % Number of UEs
OFDM.bps = 4;                   % Bits per QAM symbol (ADC/DAC resolution)
OFDM.nfft = 64;                 % FFT length
OFDM.subs = 120e3;              % Subcarrier spacing (B/M)
OFDM.BW = OFDM.nfft*OFDM.subs;  % System Bandwidth

%%% Channel parameters %%%
Chan.delay_spread = 1000e-9;    % Delay Spread of the channel
Chan.doppler = 0;
Chan.fc = 2e9;                  % center frequency
Chan.LSpeed = physconst('LightSpeed');
Chan.lambda = Chan.LSpeed/Chan.fc;
Chan.ChannelType = 'Custom';
Chan.LoS = true;
Chan.LoSKfactor = 1;
Chan.numClusters = 8;
Chan.Rays = 5;
Chan.NoiseFigure = 6; %[dB]
% Chan.NoisePowerdBm = -174 + 10*log10(OFDM.BW) + Chan.NoiseFigure;
% Chan.NoisePower = db2pow(Chan.NoisePowerdBm-30);
% Chan.NoisePower = 0.005;

%%% Exponential power delay profile
L = 40;
Tm = (1/OFDM.BW)/Chan.delay_spread;
Chan.pathDelays = (1/OFDM.BW)*([0:L-1]);
% PDP_EXP_norm = exp((-1*Tm)*([0:params.sys.L-1]));
PDP_EXP = Tm*exp((-1*Tm)*([0:L-1]));
PDP_EXP_norm = reNormalize(PDP_EXP);
Chan.pathGains=pow2db(PDP_EXP_norm);
%%%

% Chan.pathGains  = [0 -1 -2 -3 -8 -17.2 -20.8];
Chan.Noise = true; % Check false for Noise free transmission
Chan.Estim = true; % Check false for Perfect channel estimation

%%% BS and UE parameters; LoS angle calculation %%%
BS.nAntenna = Nt;               % Number of transmit antennas (Nt)
BS.AntennaType = 'ULA';         % ULA, UPA
BS.Loc = [0; 0];
UE.nAntenna = Nr;               % Number of receive antennas (Nr)
UE.AntennaType = 'ULA';         % ULA, UPA
UE.Loc = [100; 100];


