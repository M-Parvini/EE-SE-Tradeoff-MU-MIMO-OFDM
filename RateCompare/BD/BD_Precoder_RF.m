function [FRF] = BD_Precoder_RF(H, Chan, OFDM, BS, UE, WRF)
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

Hint = zeros(U*NRF,Nt);

for i=1:U
    Hk = mean(H(:,:,:,i), 3);
    wrf = WRF(:,:,i);
    Hint((i-1)*NRF+1:i*NRF,:) = wrf'*Hk;
end
coeff = 1/sqrt(Nt);
FRF = coeff*exp(1i*angle(Hint'));
end