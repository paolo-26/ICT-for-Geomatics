close all; clc; clear all;

% Satellite position.
X = 15487292.829;
Y =  6543538.932;
Z = 20727274.429;

% Parameters.
a = 6378137;  % Semi-major axis
f = 1 / 298.257223;  % Flattening

% Station coordinates.
phi = dms2degrees([45, 3, 48.114]);  % Latitude in degrees
lambda = dms2degrees([7, 39, 40.605]);  % Longitude in degrees
phi = deg2rad(phi);  % Latitude in radiants
lambda = deg2rad(lambda);  % Longitude in radiants

e2 = (2*f - f^2);
W = sqrt(1 - e2*(sin(phi))^2);

% Station position.
Xi = a * cos(phi) * cos(lambda) / W;
Yi = a * cos(phi) * sin(lambda) / W;
Zi = a * (1-e2)   * sin(phi)    / W;


R = [-sin(lambda),           cos(lambda),            0;...
     -sin(phi)*cos(lambda),  -sin(phi)*sin(lambda),  cos(phi);...
      cos(phi)*cos(lambda),  cos(phi)*sin(lambda),   sin(phi)];
M = [X-Xi; Y-Yi; Z-Zi];
enu = R * M;

% Local coordinates.
e = enu(1)
n = enu(2)
u = enu(3)
azimuth = atan(u / sqrt(n^2 + e^2));  % Azimuth in radiants
azimuth = rad2deg(azimuth)  % Azimuth in degrees.
elevation = atan(e / n);  % Azimuth in radiants
elevation = rad2deg(elevation)  % Elevation in degrees.

           