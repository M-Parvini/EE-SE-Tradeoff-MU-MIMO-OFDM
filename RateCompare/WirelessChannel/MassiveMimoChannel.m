function [Chan,H_fc,Atfc,Arfc,gains] = MassiveMimoChannel(Chan, OFDM, BS, UE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% central frequency fc
fc = Chan.fc;
% bandwidth
fs = OFDM.BW;
% OFDM subcarrier numbers
N = OFDM.nfft;
% Number of clusters
Nc = Chan.numClusters;
% Number of rays
Nray = Chan.Rays;
% BS location
BSLoc = BS.Loc;
% UE location
UELoc = UE.Loc;
d = norm(abs(UELoc-BSLoc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angle_sigma = 10/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
gamma = sqrt((Nt*Nr)/(Nc*Nray)); %normalization factor
sigma = 1; %according to the normalization condition of the H
tmax = Chan.delay_spread;

%%% Frequency-domain channel
H_fc=zeros(Nr,Nt,N);

%%% BS and UE AoD and AoA, 2 refers to angles in both azimuth and elevation
Alpha_BS = zeros(Nc*Nray, 2);
Alpha_UE = zeros(Nc*Nray, 2);
gains = zeros(Nc*Nray, 1);
delays = rand(Nc*Nray,1)*tmax;

%%% Array responses
Atfc = zeros(Nt, Nc*Nray);
Arfc = zeros(Nr, Nc*Nray);

%%% Channel modeling

for c = 1:Nc
    % Mean of the angles in azimuth and elevation
    AoD_m = unifrnd(0,pi,1,2)-pi/2;
    AoA_m = unifrnd(0,pi,1,2)-pi/2;
    
    % Azimuth and elevation angles for each cluster according to the mean
    % values from above --> according to laplacian distribution
    AoD(1,:) = laprnd(1,Nray,AoD_m(1),angle_sigma);
    AoD(2,:) = laprnd(1,Nray,AoD_m(2),angle_sigma);
    AoA(1,:) = laprnd(1,Nray,AoA_m(1),angle_sigma);
    AoA(2,:) = laprnd(1,Nray,AoA_m(2),angle_sigma);
    % AoD(1,:) = ones(1,Nray)*AoD_m(1);
    % AoD(2,:) = ones(1,Nray)*AoD_m(2);
    % AoA(1,:) = ones(1,Nray)*AoA_m(1);
    % AoA(2,:) = ones(1,Nray)*AoA_m(2);
    
    %%% Saving the angles (Used mainly in TTD part)
    Alpha_BS((c-1)*Nray+1:c*Nray,:) = AoD.';
    Alpha_UE((c-1)*Nray+1:c*Nray,:) = AoA.';

    for j = 1:Nray
        temp = (c-1)*Nray+j;

        %%% Frequency independent array responses
        Atfc(:,temp) = ArrayResponse([AoD(1,j),AoD(2,j)],BS,Chan,Chan.LSpeed/fc);
        Arfc(:,temp) = ArrayResponse([AoA(1,j),AoA(2,j)],UE,Chan,Chan.LSpeed/fc);
        %%% channel gain
        alpha = normrnd(0,sqrt(sigma/2)) + 1i*normrnd(0,sqrt(sigma/2));
        % alpha=1;
        gains(temp, 1) = alpha;
        for k=1:N
            fk = fs/(N)*(k-1-(N-1)/2);
            f = fc+fk;
            PLval = Chan.LSpeed/(4*pi*f*d);  % pathlosss model
            PLval = 1; % to cancel the impact of pathloss and focus on Small-scale fading
            H_fc(:,:,k) = H_fc(:,:,k)+alpha*PLval*...
            Arfc(:,temp)*Atfc(:,temp)'*...
            exp(-1j*2*pi*delays(temp,1)*f);
            
            % frequency flat channel
            % H_fc(:,:,k) = H_fc(:,:,k)+alpha*PLval*...
            % Arfc(:,temp)*Atfc(:,temp)';

        end
    end
end

H_fc = H_fc.*gamma;

%%% Saving the array response
% Chan.Atfc = Atfc;
% Chan.Arfc = Arfc;
Chan.normalizedPathGains = gains;

end
