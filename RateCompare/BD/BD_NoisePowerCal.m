function R = BD_NoisePowerCal(Chan, H, OFDM, FRF, FBB, WRF, WBB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of RF chains
UENRF = OFDM.UERFchain;
BSNRF = OFDM.BSRFchain;
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

R = zeros(Ns,Ns,K,U);
for k=1:K
    for j=1:U
        for i=1:U
            if j~=i
                R(:,:,k,j) = R(:,:,k,j) + ...
                        WBB(:,:,k,j)'*WRF(:,:,j)'*H(:,:,k,j)*FRF*...
                        FBB(:,:,k,i)*FBB(:,:,k,i)'*...
                        FRF'*H(:,:,k,j)'*WRF(:,:,j)*WBB(:,:,k,j);
                % R(:,:,k,j) = R(:,:,k,j) + (Pin/(U*Ns))*...
                %         WBB(:,:,k,j)'*WRF(:,:,j)'*H(:,:,k,j)*FRF*...
                %         FBB(:,:,k,i)*FBB(:,:,k,i)'*...
                %         FRF'*H(:,:,k,j)'*WRF(:,:,j)*WBB(:,:,k,j);
            end
        end
        R(:,:,k,j) = R(:,:,k,j) + sigma*WBB(:,:,k,j)'*WRF(:,:,j)'*...
                                WRF(:,:,j)*WBB(:,:,k,j);
    end
end

end