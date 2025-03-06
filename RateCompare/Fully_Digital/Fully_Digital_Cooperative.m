function SE = Fully_Digital_Cooperative(H_total, Chan, OFDM, BS, UE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notice that if the inter-user interference is zero, then the system
% simplifies to U independent streams which here we calculate one streams
% capacity --> Users have similar capacities when there is no interference
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

H = H_total(:,:,:,1);
WD = zeros(Nr,Ns,K);
FD = zeros(Nt,Ns,K);

for k=1:K
    [U_,~,V_] = svd(H(:,:,k));
    WD(:,:,k) = U_(:,1:Ns);
    FD(:,:,k) = V_(:,1:Ns);
end

%% Capacity calculation
SE=0;
% Noise power calculation [Fully connected Architecture]
Pn_Full = zeros(Ns,Ns,K);
for k=1:K
    Pn_Full(:,:,k) = sigma*WD(:,:,k)'*WD(:,:,k);
end


for k=1:K
    Pu_Full = WD(:,:,k)'*H_total(:,:,k)*...
              FD(:,:,k)*FD(:,:,k)'*...
              H_total(:,:,k)'*WD(:,:,k);

    SE = SE + ...
        log2(det(eye(Ns) + pinv(Pn_Full(:,:,k))*Pu_Full))/OFDM.nfft;
end









% H_l = zeros(U*Nr,Nt,K);
% for k=1:K
%     emptyMat = [];
%     q=num2cell(squeeze(H_total(:,:,k,:)),[1,2]);
%     for u=1:U
%         emptyMat = [emptyMat, q{u}.'];
%     end
%     H_l(:,:,k) = emptyMat.';
% end
% 
% %%% Precoders and Combiners
% WD = zeros(Nr,Ns,K,U);
% FD = zeros(Nt,Ns,K,U);
% WD_T = zeros(U*Nr, U*Ns, K);
% FD_T = zeros(Nt, U*Ns, K);
% 
% for k=1:K
%     [U_,~,V_] = svd(H_l(:,:,k));
%     WD_T(:,:,k) = U_(:,1:U*Ns);
%     FD_T(:,:,k) = V_(:,1:U*Ns);
% end
% 
% for k=1:K
%     for u=1:U
%         FD(:,:,k,u) = FD_T(:,(u-1)*Ns+1:u*Ns,k);
%     end
% end
% 
% for k=1:K
%     for u=1:U
%         WD(:,:,k,u) = WD_T((u-1)*Nr+1:u*Nr,(u-1)*Ns+1:u*Ns,k);
%     end
% end
% 
% 
% %% Capacity calculation
% SE = zeros(1, U);
% 
% % Noise power calculation [Fully connected Architecture]
% Pn_Full = Fully_Digital_NoisePowerCal(Chan, H_total, OFDM, FD, WD);
% 
% for i = 1:U
%     for k=1:K
%         Pu_Full = WD(:,:,k,i)'*H_total(:,:,k,i)*...
%                   FD(:,:,k,i)*FD(:,:,k,i)'*...
%                   H_total(:,:,k,i)'*WD(:,:,k,i);
% 
%         SE(1,i) = SE(1,i) + ...
%             log2(det(eye(Ns) + pinv(Pn_Full(:,:,k,i))*Pu_Full))/OFDM.nfft;
%     end
% end

end