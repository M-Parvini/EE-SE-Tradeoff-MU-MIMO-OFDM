function [rho_a, Ruu] = ADCModel(OFDM)
    
    b = OFDM.DAClevels;
    rho_a = QuantizationNoise(b);
    Ruu = 0;
    % RFChain = OFDM.RFchain;
    % for k=1:K
    %     Ruu(:,:,k) = Pin/Ns * FBB(:,:,k)*FBB(:,:,k)'*eye(RFChain);
    % end
        
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