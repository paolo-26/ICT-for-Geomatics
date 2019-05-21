clc; close all; clear;

% Cartographic coordinates.
east = 470139.66;
north = 5031468.37;

% WGS84 ellipsoid parameters.
e = sqrt(0.006694379990141);
e_ = sqrt(0.006739496742276);
a = 6378137;
c = 6.356752314245179e6;

% Variables.
Rp = a^2 / c;
A1 = 1 - e^2/4 - 3*e^4/64 - 5*e^6/256;

E1 = (a-c) / (a+c);

B2 = (3 * E1 / 2) - (27 * E1^3 / 32);
B4 = (21 * E1^2 / 16) - (55 * E1^4 / 32);
B6 = 151 * E1^3 / 96;
B8 = 1097 * E1^4 / 512;

mc = 0.9996;
falseEast = 5e5;
y = north / mc;
x = (east-falseEast) / mc;
lambdaMc = deg2rad(9);

theta = y / (a * A1);
xi = theta + B2*sin(2*theta) + B4*sin(4*theta) + ...
     B6*sin(6*theta) + B8*sin(8*theta);
nu =  sqrt(1 + e_^2*cos(xi)^2);
lambda1 = atan(nu * sinh(x/Rp) / cos(xi));

phi = atan(tan(xi)* cos(nu * lambda1));

% Geographic coordinates.
lambda = rad2deg(lambda1 + lambdaMc);
phi = rad2deg(phi);

array2table([phi lambda], 'VariableNames', {'phi', 'lambda'})

