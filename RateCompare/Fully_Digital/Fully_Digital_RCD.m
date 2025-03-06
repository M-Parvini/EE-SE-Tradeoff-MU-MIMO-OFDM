function SE = Fully_Digital_RCD(H_total, Chan, OFDM, BS, UE)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Precoders and Combiners
WD = zeros(Nr,Ns,K,U);
FD = zeros(Nt,Ns,K,U);

[Nr,Nt,K,U] = size(H_total);    

for m = 1:K
    for j = 1:U
        H_tildeJ = [];
        for k = setdiff(1:U,j)
            H_tildeJ = [H_tildeJ(:,:);   H_total(:,:,m,k)]; 
        end

        Fa = pinv(H_tildeJ'*H_tildeJ + ((Nr*sigma))*eye(Nt));
        [U_,~,V_] = svd(H_total(:,:,m,j)*Fa);
        
        Fb = V_(:,1:Ns);
        FD(:,:,m,j) = Fa*Fb;
        FD(:,:,m,j) = FD(:,:,m,j)/norm(FD(:,:,m,j),'fro');
        WD(:,:,m,j) = U_(:,1:Ns);
        WD(:,:,m,j) = WD(:,:,m,j)/norm(WD(:,:,m,j),'fro');
    end
end

%% Capacity calculation
SE = zeros(1, U);

% Noise power calculation [Fully connected Architecture]
Pn_Full = Fully_Digital_NoisePowerCal(Chan, H_total, OFDM, FD, WD);

for i = 1:U
    for k=1:K
        Pu_Full = WD(:,:,k,i)'*H_total(:,:,k,i)*...
                  FD(:,:,k,i)*FD(:,:,k,i)'*...
                  H_total(:,:,k,i)'*WD(:,:,k,i);
    
        SE(1,i) = SE(1,i) + ...
            log2(det(eye(Ns) + pinv(Pn_Full(:,:,k,i))*Pu_Full))/OFDM.nfft;
    end
end

end