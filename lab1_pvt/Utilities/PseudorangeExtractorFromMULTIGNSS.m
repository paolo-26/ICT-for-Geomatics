function [ RHO ] = PseudorangeExtractorFromMULTIGNSS(filePath,refPosition,error_std)
%PseudorangeExtractorFromMULTIGNSS extract pseudoranges for all available
%satellites and constellations
%   Detailed explanation goes here


settings.sigmaRho=error_std;                % [meters]
settings.observationTime=100;       % [seconds]
settings.dataSetRepository=filePath;

ECEF_USR=llh2xyz(refPosition(1),refPosition(2),refPosition(3));

load(settings.dataSetRepository);
dataSet=out(1:settings.observationTime);

%% Multi-constellation (GPS,GLONASS,BEIDOU ECEF matrix extraction
for time_index=1:settings.observationTime
SV_ECEF.GPS(:,:,time_index)=out(time_index).GPS.xyz';
SV_ECEF.GLO(:,:,time_index)=out(time_index).GLO.xyz';
SV_ECEF.BEI(:,:,time_index)=out(time_index).BEI.xyz';
SV_ECEF.GAL(:,:,time_index)=out(time_index).GAL.xyz';
end

%% Pseudorange extraction from ECEF reference [GPS]
n_sat=size(SV_ECEF.GPS,1);

for tt=1:settings.observationTime
for kk=1:n_sat
    RHO.GPS(kk,tt)=norm(SV_ECEF.GPS(kk,:,time_index)-ECEF_USR)+settings.sigmaRho*randn();
end
end

%% Pseudorange extraction from ECEF reference [GLONASS]
n_sat=size(SV_ECEF.GLO,1);

for tt=1:settings.observationTime
for kk=1:n_sat
    RHO.GLO(kk,tt)=norm(SV_ECEF.GLO(kk,:,time_index)-ECEF_USR)+settings.sigmaRho*randn();
end
end

%% Pseudorange extraction from ECEF reference [GPS]
n_sat=size(SV_ECEF.BEI,1);

for tt=1:settings.observationTime
for kk=1:n_sat
    RHO.BEI(kk,tt)=norm(SV_ECEF.BEI(kk,:,time_index)-ECEF_USR)+settings.sigmaRho*randn();
end
end

end

