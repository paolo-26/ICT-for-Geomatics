close all; clc; clear all;

phi = dms2degrees([45, 3, 48]);
phi = deg2rad(phi);
lambda = dms2degrees([7, 39, 41]);
lambda = deg2rad(lambda);
height = 0;  % Metres.

% Parameters.
a = 6378137;  % Semi-major axis
f = 1 / 298.257223;  % Flattening

e2 = (2*f - f^2);
W = sqrt(1 - e2*(sin(phi))^2);

Xi = a * cos(phi) * cos(lambda) / W;
Yi = a * cos(phi) * sin(lambda) / W;
Zi = a * (1-e2)   * sin(phi)    / W;

sat = [22504974.806   13900127.123   -2557240.727;
       -3760396.280  -17947593.853   19494169.070;
        9355256.428  -12616043.006   21189549.365;
       23959436.524    5078878.903  -10562274.680;
       10228692.060  -19322124.315   14550804.347;
       23867142.480   -3892848.382   10941892.224;
       21493427.163  -15051899.636    3348924.156;
       14198354.868   13792955.212   17579451.054;
       18493109.722    4172695.812   18776775.463;
       -8106932.299   12484531.565   22195338.169;
        8363810.808   21755378.568   13378858.106];
  
for s = 1:length(sat(:, 1))
    X = sat(s, 1);
    Y = sat(s, 2);
    Z = sat(s, 3);
    rho = sqrt((X-Xi)^2 + (Y-Yi)^2 + (Z-Zi)^2);
    D1 = (X - Xi) / rho;
    D2 = (Y - Yi) / rho;
    D3 = (Z - Zi) / rho;
    D4 = -1;
    D(s, 1:4) = [D1 D2 D3 D4];
end

Qxx = inv(D.' * D);
Qxx_star = Qxx(1:3, 1:3);

R = [-sin(lambda),           cos(lambda),            0;...
     -sin(phi)*cos(lambda),  -sin(phi)*sin(lambda),  cos(phi);...
      cos(phi)*cos(lambda),  cos(phi)*sin(lambda),   sin(phi)];

Quu = R * Qxx_star * R.';

HDOP = sqrt(Quu(1, 1) + Quu(2, 2))
PDOP = sqrt(Quu(1, 1) + Quu(2, 2) + Quu(3, 3))
GDOP = sqrt(Quu(1, 1) + Quu(2, 2) + Quu(3, 3) + Qxx(4,4))