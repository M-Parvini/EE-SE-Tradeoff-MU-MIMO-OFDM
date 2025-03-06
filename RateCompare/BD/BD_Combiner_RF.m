function [WRF] = BD_Combiner_RF(H, Chan, OFDM, BS, UE)
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
% Codebook design
D = zeros(Nr, Nr);
delta = 2*pi/Nr;
for j=1:Nr
    D(:,j) = (1/sqrt(Nr))*exp(1i*(0:Nr-1)*(j-1)*delta);
end

WRF = zeros(Nr,NRF,U);

for i=1:U
    Hk = mean(H(:,:,:,i), 3);
    temp = sum(abs(D'*Hk),2);
    [~, idx] = sort(temp, 'descend');
    WRF(:,:,i) = D(:,idx(1:NRF));
end
