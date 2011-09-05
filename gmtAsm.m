function zActuator = gmtAsm(norm)

if nargin<0
    norm = 1;
end

radialPitchOnAxis = [...
    53.03;
    35.35;
    35.34;
    35.33;
    35.32;
    35.31;
    35.29;
    35.27;
    35.25;
    35.23;
    35.20;
    35.17;
    35.14;
    35.11];

radialPitchOffAxis = [...
    52.86;
    35.49;
    35.49;
    35.48;
    35.47;
    35.46;
    35.44;
    35.43;
    35.41;
    35.39;
    35.37;
    35.34;
    35.31;
    35.28];

nActuatorRing = [...
    9;
    15;
    21;
    27;
    33;
    39;
    45;
    51;
    57;
    63;
    69;
    75;
    81;
    87];

nRing = length(nActuatorRing);
nActuator = sum(nActuatorRing);
% asmD = 1.06;
gmt = giantMagellanTelescope;
sumRadialPitchOnAxis  = 1e-3*cumsum(radialPitchOnAxis)*gmt.segmentD/1062.80e-3;
sumRadialPitchOffAxis = 1e-3*cumsum(radialPitchOffAxis)*gmt.segmentD/1053.77e-3;

% diff(sumRadialPitchOnAxis)
% diff(sumRadialPitchOffAxis)

zActuator = zeros(nActuator,7);
for kSegment=1:7
    u = 0;
    for kRing = 1:nRing
        u = (1:nActuatorRing(kRing)) + u(end);
        o        = 2*pi*(0:nActuatorRing(kRing)-1)/nActuatorRing(kRing);
        if kRing==1
            r  = ones(1,nActuatorRing(kRing))*sumRadialPitchOnAxis(kRing);
        else
            r = ones(1,nActuatorRing(kRing))*sumRadialPitchOffAxis(kRing);
        end
        zActuator(u,kSegment) = r.*exp(1i*o) + gmt.segmentCoordinate(kSegment);
    end
end
% zActuator = zActuator(:);
% zActuator(abs(zActuator)<giantMagellanTelescope.centralObscurationD/2) = [];
zActuator = zActuator/norm;

figure
h = polar(angle(zActuator),abs(zActuator),'r.');
set(h,'markersize',10)
% polar(o,rOnAxis,'rx')
% hold 
% polar(o,rOffAxis,'b+')

