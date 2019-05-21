clc; clear; close all;
format long

% Geographic coordinates.
phi = [deg2rad(dms2degrees([45 3 45.717]));
       deg2rad(dms2degrees([38 32 34.649]))];  % Radiants
lambda = [deg2rad(dms2degrees([7 47 26.292]));
          deg2rad(dms2degrees([16 50 6.493]))];  % Radiants

% Reference system: WGS84 - UTM 32 and UTM 33.
lambdaMc = deg2rad([9; 15]); % Central meridian: 9° and 15°

% WGS84 ellipsoid parameters.
e = sqrt(0.006694379990141);
e1 = sqrt(0.006739496742276);
a = 6378137;
c = 6.356752314245179e6;
mc = 0.9996;

% A.
A1 = 1 - e^2/4 - 3*e^4/64 - 5*e^6/256;
A2 = 3*e^2/8 + 3*e^4/32 + 45*e^6/1024;
A4 = 15*e^4/256 + 45*e^6/1024;
A6 = 35*e^6/3072;

% Variables.
lambda1 = (lambda - lambdaMc);
Rp = a^2 / c;
nu1 = sqrt(1+e1^2.*(cos(phi)).^2);
xi = atan(tan(phi)./cos(nu1.*lambda1));
nu =  sqrt(1+e1^2*(cos(xi)).^2);
falseEast = 5e5;

% Equations.
x = Rp * asinh(cos(xi).*tan(lambda1)./nu);
y = a*(A1*xi - A2*sin(2*xi) + A4*sin(4*xi) - A6*sin(6*xi));

% Cartographic coordinates.
east = x*mc + falseEast;
north = y*mc;

array2table([north east], 'VariableNames', {'N', 'E'},...
            'RowNames', {'Point 1', 'Point 2'})
