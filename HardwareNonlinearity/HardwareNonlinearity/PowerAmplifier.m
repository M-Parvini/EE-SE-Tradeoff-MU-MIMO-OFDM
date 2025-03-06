function [Rss, beta, alpha] = PowerAmplifier(OFDM, FRF, Ruu, rho, Nt)
%%% Power Amplifier model including Rapp, Saleh, Wiener-Hammerstein

BO = db2pow(OFDM.PABO);
K = OFDM.nfft;
PAgain = sqrt(db2pow(OFDM.PAGain));

beta = PAgain^2*(1-exp(-BO));
alpha = PAgain*(1-exp(-BO)+sqrt(BO*pi)/2*erfc(sqrt(BO)));

for k=1:K
    Rss(:,:,k) = (1-rho)*FRF*Ruu(:,:,k)*FRF'*eye(Nt);
end

end