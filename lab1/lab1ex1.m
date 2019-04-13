% Code for tasks 1-2

clc; clear; close all;
format long

% NominalUEREfixed.
P1 = load('dataset_CPT_20180328T122038');
P2 = load('dataset_HEL_20180329T160900');  % Dynamic
P3 = load('dataset_LYB_20180328T121914');
P4 = load('dataset_SHA_20180328T121804');
P5 = load('dataset_STA_20180328T121529');
P2 = load('dataset_TRN_20180328T121701');

place = '     Point';  % Label for world plot

time = length(P2.RHO.GPS(1,:));
xTicks = [0 : 400 : time];

sat = zeros(1,time);
for s = 1 : time
    sat(s) = sum(not(isnan(P2.RHO.GPS(:,s))));
end

%% Recursive Least Mean Square

for n = 1 : time  % time instants
    
    % Initialization.
    xHat = [0, 0, 0, 0]; 
    
    % Variables
    rho = P2.RHO.GPS(:,n);
    visSat = find(not(isnan(P2.RHO.GPS(:,n))));  %visible satellites indexes
    H = zeros(length(visSat),4);
    H(:,4) = 1;
                         
    for k = 1:10  % Iterations of algorithm
        
        % Step 0: rhoHat and H_matrix evaluation
        for s = 1 : length(visSat)  % Iterations on satellites
            e = P2.SAT_POS_ECEF.GPS(visSat(s)).pos(n,:) - xHat(1:3);
            rhoHat(s) = sqrt(sum(e.^2));
            H(s,1:3) = [e/rhoHat(s)];
        end  
        
        % Step 1: DeltaRhoHat  
        eRhoHat = rhoHat' - rho(visSat);
        
        % Step 2: DeltaXHat
        eXHat = (inv(H.'*H))*H.'*eRhoHat;

        % Step3: xHat
        xHat = xHat + eXHat.';
        
    end
    
    G = inv(H.'*H);
    pdop(n) = sqrt(trace(G(1:end-1, 1:end-1)));
    
    pos(n,1:4) = xHat;  %pos è in x,y,z
    clear rhoHat
    clear xHat
end

for tPos = 1 : length(pos)
    coord(tPos, 1:3) = ecef2lla(pos(tPos, 1:3));  %coord è in lon,lat,alt
end

L = length(P2.RHO.GPS(:,1))
for c = 1 : L
    r = P2.RHO.GPS(c,:);
    r = r(~(isnan(r)));
    diffR = diff(r, 2);
    sigmaUERE(c)=var(diffR);
end
R = diag(sigmaUERE);

% Save position for Google Earth
filename = '.\lab1\Helsinki.m';
writeKML_GoogleEarth(filename, coord(:, 1), coord(:, 2), coord(:, 3))

%plot(errors)
%H2=H(1:end-1,1:end-1)
%pdop = sqrt(trace(inv(H2'*H2)))



%% Plot: Visible satellites vs time.
figure(1)
plot(sat, 'linewidth', 2)
xlabel("Time")
ylabel("Number of visible satellites")
title("Number of visible satellites vs Time")
ylim([min(sat)-1, max(sat)+1])
xlim([0,time])
xticks(xTicks);
grid on

%% Plot: Pseudoranges vs time.
figure(2)
plot([1:1:time], P2.RHO.GPS', '.')
xlabel("Time")
ylabel("Pseudorange")
title("Pseudoranges vs Time")
xlim([0,time])
xticks(xTicks);
grid on

%% Plot: World Map Plot
figure(3)
WorldMapPlotter(coord(time,1), coord(time,2), place, 'World')

%% Plot: Scatter plot of positions
% Real position is the mean of all positions.
realPos = [mean(pos(:,1)), mean(pos(:,2)), mean(pos(:,3))];
realPosLla = ecef2lla(realPos);
errors = realPos - pos(:,1:3);
sigmaX = sqrt(trace(cov(errors)));
fprintf("sigmaX = %.6f\n", sigmaX)
fprintf("Avg. PDOP = %.2d metres\n", mean(pdop))
fprintf("Avg. sigmaUERE = %.2d\n", nanmean(sigmaUERE))
res = mean(pdop)*(nanmean(sigmaUERE));
fprintf("PDOP * sigmaUERE = %.6f\n", res)

figure(4)
geoshow(coord(:,1), coord(:,2),'DisplayType', 'point', 'Marker', '.')
geoshow(realPosLla(1), realPosLla(2),'DisplayType', 'point',...
        'MarkerEdgeColor','blue', 'Marker','o', 'MarkerFaceColor', 'green')
legend('Estimated positions', 'Real position')
xlabel('Longitude')
ylabel('Latitude')
title('Positions')
axis square; grid on;

%% Plot: PDOP vs time.
figure(5)
hold on
plot(pdop)
title("Position dilution of precision")
grid minor
xlabel('Time')
ylabel('PDOP')
xlim([0,time])
xticks(xTicks);



std(errors)











%% WLMS

for n = 1 : time  % time instants
    
    % Initialization.
    xHat = [0, 0, 0, 0]; 
    
    % Variables
    rho = P2.RHO.GPS(:,n);
    visSat = find(not(isnan(P2.RHO.GPS(:,n))));  % Visible sat. indexes
    H = zeros(length(visSat),4);
    H(:,4) = 1;
                         
    for k = 1:10  % Iterations of algorithm
        
        % Step 0: rhoHat and H_matrix evaluation
        for s = 1 : length(visSat)  % Iterations on satellites
            e = P2.SAT_POS_ECEF.GPS(visSat(s)).pos(n,:) - xHat(1:3);
            rhoHat(s) = sqrt(sum(e.^2));
            H(s,1:3) = [e/rhoHat(s)];
        end  
        
        % Step 1: DeltaRhoHat  
        eRhoHat = rhoHat' - rho(visSat);
  
        % Step 2: H_w
        W = inv(R(visSat, visSat));
        H_w = inv(H.'*W*H)*H.'*W;
        % Step 3: DeltaXHat
        eXHat = H_w * eRhoHat; 
        
        % Step 4: xHat
        xHat = xHat + eXHat.';
                
        
    end
    
    H_w = H_w.';
    G = inv(H_w.'*H_w);
    pdop(n) = sqrt(trace(G(1:end-1, 1:end-1)));
    
    pos(n,1:4) = xHat;  %pos è in x,y,z
    clear rhoHat
    clear xHat
end



%% Plot: Visible satellites vs time.
figure(11)
plot(sat, 'linewidth', 2)
xlabel("Time")
ylabel("Number of visible satellites")
title("Number of visible satellites vs Time")
ylim([min(sat)-1, max(sat)+1])
xlim([0,time])
xticks(xTicks);
grid on

%% Plot: Pseudoranges vs time.
figure(12)
plot([1:1:time], P2.RHO.GPS', '.')
xlabel("Time")
ylabel("Pseudorange")
title("Pseudoranges vs Time")
xlim([0,time])
xticks(xTicks);
grid on

%% Plot: World Map Plot
figure(13)
WorldMapPlotter(coord(time,1), coord(time,2), place, 'World')

%% Plot: Scatter plot of positions
% Real position is the mean of all positions.
realPos = [mean(pos(:,1)), mean(pos(:,2)), mean(pos(:,3))];
realPosLla = ecef2lla(realPos);
errors = realPos - pos(:,1:3);
sigmaX = sqrt(trace(cov(errors)));

fprintf("Avg. PDOP = %.2d metres\n", mean(pdop))
fprintf("Error: %.2d metres\n", sigmaX)
fprintf("Avg. sigmaUERE = %.2d\n", nanmean(sigmaUERE))
mean(pdop)*(nanmean(sigmaUERE))

figure(14)
geoshow(coord(:,1), coord(:,2),'DisplayType', 'point', 'Marker', '.')
geoshow(realPosLla(1), realPosLla(2),'DisplayType', 'point',...
        'MarkerEdgeColor','blue', 'Marker','o', 'MarkerFaceColor', 'green')
legend('Estimated positions', 'Real position')
xlabel('Longitude')
ylabel('Latitude')
title('Positions')
axis square; grid on;

%% Plot: PDOP vs time.
figure(5)
plot(pdop)
title("Position dilution of precision")
legend
xlabel('Time')
ylabel('PDOP')
xlim([0,time])
xticks(xTicks);



%% end
close figure 3 13
std(errors)