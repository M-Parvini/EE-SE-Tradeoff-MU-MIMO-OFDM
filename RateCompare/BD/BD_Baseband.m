function [FBB, WBB] = BD_Baseband(H, Chan, OFDM, BS, UE, FRF, WRF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% Number of RF chains
UE_rf = OFDM.UERFchain;
BS_rf = OFDM.BSRFchain;
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
% baseband precoder and combiner
FBB = zeros(BS_rf,Ns,K,U);
FBB_T = [];
WBB = zeros(UE_rf,Ns,K,U);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Heff = zeros(UE_rf,BS_rf,K,U);
Vnull = zeros(BS_rf,UE_rf,K,U);
VNs = zeros(UE_rf,Ns,K,U);

for u=1:U
    for k=1:K
        Wk = WRF(:,:,u);
        Hk = H(:,:,k,u);
        Heff(:,:,k,u) = Wk'*Hk*FRF;
    end
end


%%%
for k=1:K
    for u=1:U
        Hbar = [];
        for j=1:U
            if j~=u
                Heq_tilde = Heff(:,:,k,j);
                Hbar = [Hbar,Heq_tilde.'];
            end
        end
        [~,St,Vtilde] = svd(Hbar.');
        rt = nnz(St);
        Vnull(:,:,k,u) = Vtilde(:,rt+1:end);
    end
end

%%%

% 
% 
% for k=1:K
%     for u=1:U
%         Hbar = [];
%         for j=1:U
%             if j~=u
%                 Heq_tilde = Heff(:,:,k,j);
%                 Hbar = [Hbar,Heq_tilde.'];
%             end
%         end
%         Vnull(:,:,k,u) = null(Hbar.');
%     end
% end

for k=1:K
    for u=1:U
        [S,~,D] = svd(Heff(:,:,k,u)*Vnull(:,:,k,u));
        WBB(:,:,k,u) = S(:,1:Ns);
        VNs(:,:,k,u) = D(:,1:Ns);
    end
end


% Normalization
for k=1:K
    for u=1:U
        FBB(:,:,k,u) = Vnull(:,:,k,u)*VNs(:,:,k,u);
        FBB(:,:,k,u) = FBB(:,:,k,u)/(norm(FRF*FBB(:,:,k,u),'fro'));
    end
    for u = 1:U
        WBB(:,:,k,u) = WBB(:,:,k,u) / norm(WRF(:,:,u) * WBB(:,:,k,u),'fro');
    end
end

end





% % Normalization
% for k=1:K
%     FBB_T = [];
%     for i=1:U
%         FBB(:,:,k,i) = Vnull(:,:,k,i)*VNs(:,:,k,i);
%         FBB_T = [FBB_T, FBB(:,:,k,i)];
%     end
%     FBB_T = sqrt(U*Ns) * FBB_T / norm(FRF * FBB_T,'fro');
%     for i = 1:U
%         WBB(:,:,k,i) = sqrt(Ns) * WBB(:,:,k,i) / norm(WRF(:,:,i) * WBB(:,:,k,i),'fro');
%         FBB(:,:,i) = FBB_T(:,(i-1)*Ns+1:i*Ns);
%     end
% end