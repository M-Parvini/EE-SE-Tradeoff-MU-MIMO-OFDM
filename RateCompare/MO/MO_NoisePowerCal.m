function R = MO_NoisePowerCal(Chan, H, OFDM, FRF, FBB, WD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of RF chains
UENRF = OFDM.UERFchain;
BSNRF = OFDM.BSRFchain;
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

R = zeros(Ns,Ns,K,U);
for k=1:K
    for j=1:U
        for i=1:U
            if j~=i
                R(:,:,k,j) = R(:,:,k,j) + ...
                        WD(:,:,k,j)'*H(:,:,k,j)*FRF*...
                        FBB(:,:,k,i)*FBB(:,:,k,i)'*...
                        FRF'*H(:,:,k,j)'*WD(:,:,k,j);
                % R(:,:,k,j) = R(:,:,k,j) + (Pin/(U*Ns))*...
                %         WBB(:,:,k,j)'*WRF(:,:,j)'*H(:,:,k,j)*FRF*...
                %         FBB(:,:,k,i)*FBB(:,:,k,i)'*...
                %         FRF'*H(:,:,k,j)'*WRF(:,:,j)*WBB(:,:,k,j);
            end
        end
        R(:,:,k,j) = R(:,:,k,j) + sigma*WD(:,:,k,j)'*WD(:,:,k,j);
    end
end

end