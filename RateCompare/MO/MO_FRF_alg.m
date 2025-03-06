function FRF = MO_FRF_alg(FRF_ini, Heff_l, phi_l, W_l, T_l, Chan, OFDM, BS, UE)
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
% H_total = zeros(Nr,Nt,K,U);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the manifold space of the analog precoder
frf_manifold = complexcirclefactory(Nt*Nt_rf);
problem.M = frf_manifold;


problem.cost = @(x)FRF_Cost(x,Heff_l, phi_l, W_l, T_l, Nt, Nt_rf,K);
problem.egrad = @(x)FRF_Grad(x,Heff_l, phi_l, W_l, T_l, Nt, Nt_rf,K);

[x,iter] = conjugategradient(problem,FRF_ini(:));
FRF = reshape(x,Nt,Nt_rf);


end