function [SE] = MO_Hybrid(H_total, Chan, OFDM, BS, UE)
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
InitializeHelp = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization of precoders and combiners
FRF = zeros(Nt,Nt_rf);
FBB = zeros(Nt_rf,Ns,K,U);
WRF = zeros(Nr,Nr_rf,U);
WBB = zeros(Nr_rf,Ns,K,U);
WD  = zeros(Nr,Ns,K,U); % W_D = WRF*WBB
WD_ini  = zeros(Nr,Ns,K,U);
%%% Initialization of weights for MMSE problem
T = zeros(Ns,Ns,K,U); % (\Gamma in latex and paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% algorithm phase %%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization phase
% Initialization for FRF
for i = 1 : Nt_rf %initialize V_RF
    FRF_ini(:,i) = exp(1i*unifrnd(0,2*pi,Nt,1));
end
if InitializeHelp
    FRF_ini = BD_FRF;
end

% Initialization for WD
for u = 1 : U
    for k = 1 : K
        if InitializeHelp
            WD_ini(:,:,k,u) = BD_WRF(:,:,u)*BD_WBB(:,:,k,u);
        else
            WD_ini(:,:,k,u) = exp(1i*unifrnd(0,2*pi,Nr,Ns));
        end
    end
end
% Initialization for T (\Gamma in latex and paper)
for u = 1 : U
    for k = 1 : K
        T(:,:,k,u) = eye(Ns);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining and initializing the variables
W_l = zeros(U*Nr,U*Ns,K);
for k=1:K
    q=num2cell(squeeze(WD_ini(:,:,k,:)),[1,2]);
    W_l(:,:,k) = blkdiag(q{:});
end

H_l = zeros(U*Nr,Nt,K);
for k=1:K
    emptyMat = [];
    q=num2cell(squeeze(H_total(:,:,k,:)),[1,2]);
    for u=1:U
        emptyMat = [emptyMat, q{u}.'];
    end
    H_l(:,:,k) = emptyMat.';
end

T_l = zeros(U*Ns,U*Ns,K);
for k=1:K
    q=num2cell(squeeze(T(:,:,k,:)),[1,2]);
    T_l(:,:,k) = blkdiag(q{:});
end

phi_l = zeros(1,K);
beta_l = zeros(1,K);
for k=1:K
    phi_l(1,k) = sigma*trace(T_l(:,:,k)*W_l(:,:,k)'*W_l(:,:,k));
end

Heff_l = zeros(U*Ns,Nt,K);
for k=1:K
    Heff_l(:,:,k) = W_l(:,:,k)'*H_l(:,:,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Iterative algorithm initiate
iter = 0;
Delta = 1;
new_wmmse = U*Ns+1;

while (iter<15 && Delta>=1e-4)
    % Compute FRF based on MO method
    [FRF] = MO_FRF_alg(FRF_ini, Heff_l, phi_l, W_l, T_l, Chan, OFDM, BS, UE);
    FRF_ini = FRF;
    % Compute beta and FBB
    FBB_bar = zeros(Nt_rf,Nt_rf,K);
    for k=1:K
        FBB_bar(:,:,k) = FRF'*Heff_l(:,:,k)'*T_l(:,:,k)*Heff_l(:,:,k)*FRF+...
                         phi_l(1,k)*FRF'*FRF;
        beta_l(1,k) = 1/(norm(FRF*FBB_bar(:,:,k)^(-1)*FRF'*Heff_l(:,:,k)'*T_l(:,:,k),'fro'));
        FBB_D(:,:,k) = beta_l(1,k)*FBB_bar(:,:,k)^(-1)*FRF'*Heff_l(:,:,k)'*T_l(:,:,k);
        % FBB = zeros(Nt_rf,Ns,K,U);
        for u=1:U
            FBB(:,:,k,u) = FBB_D(:,(u-1)*Ns+1:u*Ns,k);
        end
    end

    % Compute WD
    % interference
    Interference = zeros(Nr,Nr,K,U);
    for k=1:K
        for j=1:U
            for i=1:U
                if i~=j
                    Interference(:,:,k,j)=Interference(:,:,k,j)+...
                        H_total(:,:,k,j)*FRF*FBB(:,:,k,i)*FBB(:,:,k,i)'*FRF'*...
                        H_total(:,:,k,j)';
                end
            end
            Interference(:,:,k,j) = Interference(:,:,k,j) + sigma*eye(Nr);
        end
    end
    % WD
    for k=1:K
        for u=1:U
            WD(:,:,k,u) = beta_l(1,k)*(Interference(:,:,k,u)+...
                H_total(:,:,k,u)*FRF*FBB(:,:,k,u)*FBB(:,:,k,u)'*FRF'*...
                H_total(:,:,k,u)')^(-1)*H_total(:,:,k,u)*FRF*FBB(:,:,k,u);
        end
    end
    % T
    for k=1:K
        for u=1:U
            T(:,:,k,u) = eye(Ns)+FBB(:,:,k,u)'*FRF'*H_total(:,:,k,u)'*...
                Interference(:,:,k,u)^(-1)*H_total(:,:,k,u)*FRF*FBB(:,:,k,u);
        end
    end
    % for k=1:K
    %     for u=1:U
    %         T(:,:,k,u) = ((eye(Ns)-(1/beta_l(1,k))*WD(:,:,k,u)'*H_total(:,:,k,u)*...
    %             FRF*FBB(:,:,k,u)-(1/beta_l(1,k))*FBB(:,:,k,u)'*FRF'*...
    %             H_total(:,:,k,u)'*WD(:,:,k,u)+(1/(beta_l(1,k)^2))*WD(:,:,k,u)'*...
    %             H_total(:,:,k,u)*FRF*FBB(:,:,k,u)*FBB(:,:,k,u)'*FRF'*...
    %             H_total(:,:,k,u)'*WD(:,:,k,u))* + ...
    %             (1/(beta_l(1,k)^2))*WD(:,:,k,u)'*Interference(:,:,k,u)*...
    %             WD(:,:,k,u)+(1/(beta_l(1,k)^2))*sigma*WD(:,:,k,u)'*WD(:,:,k,u))^(-1);
    %     end
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-initialization
    W_l = zeros(U*Nr,U*Ns,K);
    for k=1:K
        q=num2cell(squeeze(WD(:,:,k,:)),[1,2]);
        W_l(:,:,k) = blkdiag(q{:});
    end
    
    T_l = zeros(U*Ns,U*Ns,K);
    for k=1:K
        q=num2cell(squeeze(T(:,:,k,:)),[1,2]);
        T_l(:,:,k) = blkdiag(q{:});
    end
    
    phi_l = zeros(1,K);
    for k=1:K
        phi_l(1,k) = sigma*trace(T_l(:,:,k)*W_l(:,:,k)'*W_l(:,:,k));
    end
    
    Heff_l = zeros(U*Ns,Nt,K);
    for k=1:K
        Heff_l(:,:,k) = W_l(:,:,k)'*H_l(:,:,k);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% MMSE Calculation
    for i=1:K
        MMSE(i)=trace(eye(U*Ns))-log(det(T_l(:,:,i)));
    end
    old_wmmse = new_wmmse;
    new_wmmse = mean(MMSE);
    Delta = old_wmmse - new_wmmse;
    iter = iter+1;

end


%% Capacity calculation
SE = zeros(1, U);

% Noise power calculation [Fully connected Architecture]
Pn_Full = MO_NoisePowerCal(Chan, H_total, OFDM, FRF, FBB, WD);

for i = 1:U
    for k=1:K
        % Pu_Full = (OFDM.Pin/(U*Ns))*...
        %                 WBB(:,:,k,i)'*WRF(:,:,i)'*H_total(:,:,k,i)*FRF*...
        %                 FBB(:,:,k,i)*FBB(:,:,k,i)'*...
        %                 FRF'*H_total(:,:,k,i)'*WRF(:,:,i)*WBB(:,:,k,i);

        Pu_Full = WD(:,:,k,i)'*H_total(:,:,k,i)*FRF*...
                  FBB(:,:,k,i)*FBB(:,:,k,i)'*...
                  FRF'*H_total(:,:,k,i)'*WD(:,:,k,i);
    
        SE(1,i) = SE(1,i) + ...
            log2(det(eye(Ns) + pinv(Pn_Full(:,:,k,i))*Pu_Full))/OFDM.nfft;
    end
end


end