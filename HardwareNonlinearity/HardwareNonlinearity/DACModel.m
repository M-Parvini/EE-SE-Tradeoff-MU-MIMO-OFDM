function [rho_d, Ruu] = DACModel(OFDM, FBB)
    U = OFDM.nUEs;
    Pin = OFDM.Pin;
    K  = OFDM.nfft;
    Ns = OFDM.numStreams;
    b = OFDM.DAClevels;
    rho_d = QuantizationNoise(b);
    RFChain = OFDM.BSRFchain;
    Ruu = zeros(RFChain,RFChain,K);
    for k=1:K
        for u=1:U
            Ruu(:,:,k) = Ruu(:,:,k) + Pin/(U*Ns) *...
                         FBB(:,:,k,u)*FBB(:,:,k,u)'*eye(RFChain);
            % Ruu(:,:,k) = Ruu(:,:,k) + Pin/(U*Ns) * ...
            %              FBB(:,:,k,u)*FBB(:,:,k,u)'*eye(RFChain);
        end
    end


        
end

function rho = QuantizationNoise(b)
    switch b
        case 1
            rho = 0.3634;
        case 2
            rho = 0.1175;
        case 3
            rho = 0.03454;
        case 4
            rho = 0.009497;
        case 5
            rho = 0.002499;
        otherwise
            rho = 0.5*pi*sqrt(3)*2^(-2*b);
    end

end