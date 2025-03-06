function R = LimFeed_NoisePowerCal(Chan, H, OFDM, FRF, FBB, WRF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of OFDM symbols
Ns = OFDM.numStreams;
% Number of UEs
U = OFDM.nUEs;
% % Pin
% Pin = OFDM.Pin;
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

R = zeros(K,U);
for k=1:K
    for j=1:U
        for i=1:U
            if j~=i
                R(k,j) = R(k,j) + ...
                        WRF(:,j)'*H(:,:,k,j)*FRF*...
                        FBB(:,i,k)*FBB(:,i,k)'*...
                        FRF'*H(:,:,k,j)'*WRF(:,j);
            end
        end
        R(k,j) = R(k,j) + sigma;
    end
end

end