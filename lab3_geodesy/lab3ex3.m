close all;
clear all;
clc;
format long

X = 4499525.4271;
Y = 585034.1293;
Z = 4467910.3596;
threshold = 10e-8;
e2 = 0.00672267002233323;  % Hayford
a = 6378388;  %Hayford
lambda = atan(Y/X);

%1
r = sqrt(X^2 + Y^2);
phi = atan(Z/r);

for k = 1 : 50
    
%2
W = sqrt(1 - (e2)*(sin(phi))^2);
N = a / W;

%3
h = (X / (cos(phi) * cos(lambda))) - N;
phi = atan(Z / (r * (1- ((e2*N)/(N+h)))));

% To save data.
% h_vect(k) = h;
% phi_vect(k) = phi;

end


phi
lambda
h




    
  