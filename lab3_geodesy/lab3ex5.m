close all;
clear;
clc;
format long;

% Alessandria.
Xa = 4472675.1340 / 1e6;
Ya = 677732.7296 / 1e6;
Za = 4481406.6192 / 1e6;

Xb = 4472674.8193 / 1e6;
Yb = 677733.0708 / 1e6;
Zb = 4481406.8974 / 1e6;

% Mondovì-
% Xa = 4523227.3918  / 1e6;
% Ya = 621792.1742   / 1e6;
% Za = 4439482.9654  / 1e6;
% 
% Xb = 4523227.0920  / 1e6;
% Yb = 621792.5562   / 1e6;
% Zb = 4439483.2501  / 1e6;

%TORINO-
% Xa = 4472544.5785  / 1e6;
% Ya = 601634.0764   / 1e6;
% Za = 4492545.0107  / 1e6;
% 
% Xb = 4472544.2920  / 1e6;
% Yb = 601634.4325   / 1e6;
% Zb = 4492545.2622  / 1e6;

A = [1 0 0  Xa   0    -Za    Ya;
     0 1 0  Ya  -Za   0     -Xa;  
     0 0 1  Za   Ya   Xa    0  ];
 
lo = [Xb-Xa; Yb-Ya; Zb-Za];
 
x = inv(A.'*A) * (A.'*lo);
v = (A*x - lo) *1e6;

Names = {'Tx (translation)', 'Ty (translation)', 'Tz (translation)', 'lambda (scaling)', 'Rx (rotation)', 'Ry (rotation)', 'Rz (rotation)'};
vNames = {'v1', 'v2', 'v3'};
x = array2table(x, 'rowNames', Names)
v = array2table(v, 'rowNames', vNames)