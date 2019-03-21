close all; clc; clear all;

% Satellite position
X = 15487292.829;
Y = 6543538.932;
Z = 20727274.429;

% Parameters.
a = 6378137.0;  % Semi-major axis
f = 1 / 298.257223;  % Flattening

% Station coordinates.
phi = 45.06336;
lambda = 7.661279;
phi = deg2rad(phi);
lambda = deg2rad(lambda);

e2 = (2*f - f^2);
W = sqrt(1 - e2*(sin(phi))^2);

% Station position.
Xi = a * cos(phi) * cos(lambda) / W;
Yi = a * cos(phi) * sin(lambda) / W;
Zi = a * (1-e2) * sin(phi) / W;


R = [-sin(lambda), cos(lambda), 0;...
     -sin(phi)*cos(lambda), -sin(phi)*sin(lambda), cos(phi);...
     cos(phi)*cos(lambda), cos(phi)*sin(lambda), sin(phi)];
M = [X-Xi; Y-Yi; Z-Zi];
res = R * M;

% Local coordinates.
e = res(1)
n = res(2)
u = res(3)
 
azimuth = rad2deg(atan(u / sqrt(n^2 + e^2)))
elevation = rad2deg(atan(e / n))

           