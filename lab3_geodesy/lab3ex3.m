% The algorithm is iterative and the transformation introduces some errors.
close all;
clear all;
clc;
format long

X = [4499525.4271; 4495694.2695; 4503484.7172; 4498329.3715];
Y = [585034.1293; 592457.8605; 578160.7507; 562840.7651];
Z = [4467910.3596; 4470744.7781; 4465024.3002; 4472537.6125];

threshold = 10e-8;
e2 = 0.00672267002233323;  % Hayford
a = 6378388;  % Hayford
lambda = atan(Y./X);

% Step 1.
r = sqrt(X.^2 + Y.^2);
phi = atan(Z./r);

for k = 1 : 10  % 7 iterations are needed
    % Step 2.
    W = sqrt(1 - (e2).*(sin(phi)).^2);
    N = a ./ W;

    % Step 3.
    h = (X ./ (cos(phi) .* cos(lambda))) - N;
    phi = atan(Z ./ (r .* (1 - ((e2*N)./(N+h)))));

    % To save data.
    h_mat(k, 1:4) = h;
    phi_mat(k, 1:4) = phi;
end

figure(1)
plot(h_mat)
title('Hayford')
grid minor
xlabel('Number of iterations')
ylabel('h [m]')
legend(['Point 1'; 'Point 2'; 'Point 3'; 'Point 4'], 'location', 'best')

colNames = {'phi', 'lambda', 'h'};
rowNames = {'Point1', 'Point2', 'Point3', 'Point4'};
hayford = array2table([rad2deg(phi) rad2deg(lambda) h],...
                    'VariableNames', colNames,...
                    'rowNames', rowNames)
colName = {'Point1', 'Point2', 'Point3', 'Point4'};
hHayford = array2table(h_mat, 'VariableName', colName)

e2 = 0.006694379990141;  % WGS84
a = 6378137;  % WGS84
lambda = atan(Y./X);

% Step 1.
r = sqrt(X.^2 + Y.^2);
phi = atan(Z./r);

for k = 1 : 10  % 7 iterations are needed
    % Step 2.
    W = sqrt(1 - (e2).*(sin(phi)).^2);
    N = a ./ W;

    % Step 3.
    h = (X ./ (cos(phi) .* cos(lambda))) - N;
    phi = atan(Z ./ (r .* (1 - ((e2*N)./(N+h)))));

    % To save data.
    h_mat(k, 1:4) = h;
    phi_mat(k, 1:4) = phi;
end

figure(2)
plot(h_mat)
title('WGS84')
grid minor
xlabel('Number of iterations')
ylabel('h [m]')
legend(['Point 1'; 'Point 2'; 'Point 3'; 'Point 4'], 'location', 'best')

colNames = {'phi', 'lambda', 'h'};
rowNames = {'Point1', 'Point2', 'Point3', 'Point4'};
WGS84 = array2table([rad2deg(phi) rad2deg(lambda) h],...
                    'VariableNames', colNames,...
                    'rowNames', rowNames)
colName = {'Point1', 'Point2', 'Point3', 'Point4'};
hWGS84 = array2table(h_mat, 'VariableName', colName)