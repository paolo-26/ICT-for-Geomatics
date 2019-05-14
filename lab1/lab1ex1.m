% Code for tasks 1-2

clc; clear; close all;
format short

% NominalUEREfixed.
P1n = load('dataset_CPT_20180328T122038');
P2n = load('dataset_HEL_20180328T122158');  
P3n = load('dataset_LYB_20180328T121914');
P4n = load('dataset_SHA_20180328T121804');
P5n = load('dataset_STA_20180328T121529');
P6n = load('dataset_TRN_20180328T121701');

% RealisticUEREdynamic.
P1r = load('dataset_CPT_20180329T160947');
P2r = load('dataset_HEL_20180329T160900');
P3r = load('dataset_LYB_20180329T161023');
P4r = load('dataset_SHA_20180329T161103');
P5r = load('dataset_STA_20180329T161418');
P6r = load('dataset_TRN_20180329T161139');

% Realistic 2019
% test1 = load('dataset_1_20180329T160947');
% test2 = load('dataset_2_20180329T160900');
% test3 = load('dataset_3_20180329T161023');
% test4 = load('dataset_4_20180329T161103');
% test5 = load('dataset_5_20180329T161418');
% test6 = load('dataset_CPT_20180329T160947');

place = '     Point';  % Label for world plot

P = P1n;

cnt = 0;
rhoMin = NaN;
rhoMax = NaN;
minCor = [NaN NaN];
maxCor = [NaN NaN];
for set =  ["GPS", "GAL", "GLO", "BEI"]
% rhoMin = min(rhoMin, min(min(P.RHO.(set))));
% rhoMax = max(rhoMax, max(max(P.RHO.(set))));
cnt = cnt + 1;
    
time = length(P.RHO.(set)(1,:));
xTicks = 0 : 900 : time;

sat = zeros(1,time);
for s = 1 : time
    sat(s) = sum(not(isnan(P.RHO.(set)(:,s))));
end

%% Recursive Least Mean Square

L = length(P.RHO.(set)(:,1));
sigmaUERE = zeros(1, L);
for c = 1 : L
    r = P.RHO.(set)(c,:);
    r = r(~(isnan(r)));
    diffR = diff(r, 2);
    sigmaUERE(c)=std(diffR);
end
R = diag(sigmaUERE).^2;

for n = 1 : time  % time instants
    
    % Initialization.
    xHat = [0, 0, 0, 0]; 
    
    % Variables
    rho = P.RHO.(set)(:,n);
    visSat = find(not(isnan(P.RHO.(set)(:,n))));  %visible satellites indexes
    H = zeros(length(visSat),4);
    H(:,4) = 1;
                         
    for k = 1:7  % Iterations of algorithm
        
        % Step 0: rhoHat and H_matrix evaluation
        for s = 1 : length(visSat)  % Iterations on satellites
            e = P.SAT_POS_ECEF.(set)(visSat(s)).pos(n,:) - xHat(1:3);
            rhoHat(s) = sqrt(sum(e.^2));
            H(s,1:3) = e/rhoHat(s);
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
    gdop(n) = sqrt(trace(G));
    deltaX(n) = ((inv(H.'*H))*H.')*nanmean(sigmaUERE);
    sigmaX(n) = sqrt(trace(inv(H.'*H) * H.'*R(visSat, visSat)*H*inv(H.'*H)));
    
    pos(n,1:4) = xHat;  %pos è in x,y,z
    clear rhoHat
    clear xHat
end

% posError = std(pos);
% fprintf("LMS sigmaX for %s =\n", set)
% display(posError')
fprintf("LMS sigmaX for %s = %.4f \n", set, sqrt(trace(cov(pos(1:3, 1:3)))))

cov_pos = cov(pos);
lmsSigmaXfull(cnt) = sqrt(trace(cov_pos));
cov_pos = cov_pos(1:3, 1:3);
cov_pos = trace(cov_pos);
lmsSigmaX(cnt) = sqrt(cov_pos);

for tPos = 1 : length(pos)
    coord(tPos, 1:3) = ecef2lla(pos(tPos, 1:3));  %coord è in lon,lat,alt
end

mPdop(cnt,1) = mean(pdop);

% Save position for Google Earth
% filename = '.\lab1\realistic_gps_capetown_lms';
% writeKML_GoogleEarth(filename, coord(:, 1), coord(:, 2), coord(:, 3))

%plot(errors)
%H2=H(1:end-1,1:end-1)
%pdop = sqrt(trace(inv(H2'*H2)))


%% Plot: Visible satellites vs time.

k = 1;
satPlot = zeros(6, time);
for PS = [P1n, P2n, P3n, P4n, P5n, P6n]
%for PS = [test1, test2, test3, test4, test5, test6]
%for PS = [testerer]
%for PS = [P1r, P2r, P3r, P4r, P5r, P6r]
    
    for s = 1 : time
        satPlot(k,s) = sum(not(isnan(PS.RHO.(set)(:,s))));
    end
    k = k+1;
end
figure(1)
subplot(2, 2, cnt)
hold on
plot(satPlot', 'linewidth', 1.5)
xlabel("Time")
ylabel("No. of visible satellites")
title(set)%, "No. of visible satellites over time"])
%ylim([5 20])
if set == 'GPS'
    ylim([0 24])
end
% ylim([min(sat)-1, max(sat)+1])
xlim([0,time])
xticks(xTicks);
grid on
grid minor

%% Plot: Pseudoranges vs time.
figure(2)
subplot(2, 2, cnt)
plot(P.RHO.(set)', 'linewidth', 1)
xlabel("Time")
ylabel("Pseudorange [m]")
title(set)
xlim([0,time])
xticks(xTicks);
grid on
grid minor

%% Plot: World Map Plot
% figure(3)
% hold on
% WorldMapPlotter(coord(time,1), coord(time,2), place, 'World')


% figure(99)
% hold on
% lat = coord(time,1);
% lon = coord(time,2);
% name = '      TRN';
% plotm(lat,lon,'ro','MarkerSize',10)
% plotm(lat,lon,'r+','MarkerSize',30,'LineWidth',1)
% textm(lat,lon,name)

%% Plot: Scatter plot of positions
% Real position is the mean of all positions.
realPos = [mean(pos(:,1)), mean(pos(:,2)), mean(pos(:,3))];
%realPos = mean(pos);
realPosLla = ecef2lla(realPos);
errors = realPos - pos(:,1:3);
%errors = realPos - pos;
%posLla = ecef2lla(pos(:,1:3))
%errors = realPosLla - posLla(:,1:3);
%sigmaX = sqrt(trace(cov(errors)));
% fprintf("sigmaX = %.6f\n", sigmaX);
fprintf("Avg. PDOP = %.6f metres\n", mean(pdop))
fprintf("Avg. GDOP = %.6f metres\n", mean(gdop))
fprintf("Avg. sigmaUERE = %.6f\n", nanmean(sigmaUERE))
res = mean(gdop)*(nanmean(sigmaUERE));
fprintf("PDOP * sigmaUERE = %.6f\n", res);

figure(4)
subplot(2, 4, cnt)
hold on
minCor(1) = min(minCor(1), min(coord(:,1)));
minCor(2) = min(minCor(2), min(coord(:,2)));
maxCor(1) = max(maxCor(1), max(coord(:,1)));
maxCor(2) = max(maxCor(2), max(coord(:,2)));
geoshow(coord(:,1), coord(:,2),'DisplayType', 'point', 'Marker', '.',...
        'MarkerEdgeColor', 	[0.6350, 0.0780, 0.1840])
geoshow(realPosLla(1), realPosLla(2),'DisplayType', 'point',...
        'MarkerEdgeColor','black', 'Marker','o',...
        'MarkerFaceColor', [0.4660, 0.6740, 0.1880])
legend('Estimated', 'Real')
xlabel('Longitude')
ylabel('Latitude')
%xlim([24.936 24.94])
%ylim([60.1696 60.17])
title(set)%'Positions')
axis equal; axis square; grid on;


% %% Plot: PDOP vs time.
% figure(5)
% hold on
% plot(pdop)
% title("Position dilution of precision")
% grid minor
% xlabel('Time')
% ylabel('PDOP')
% xlim([0,time])
% xticks(xTicks);

%% Plot: PDOP vs time.
figure(5)
subplot(2, 2, cnt)
hold on
% plot(gdop)
plot(gdop*nanmean(sigmaUERE))
plot(sigmaX)
title(set) %"Dilution of precision")
grid minor
xlabel('Time')
%ylabel('PDOP')
xlim([0,time])
xticks(xTicks);
legend('GDOP \cdot \sigma_{UERE}', '\sigma_X', 'location', 'best')

std(errors);
clear R  sigmaUERE sigmaX H G
end

% 
% figure(2)
% for k = 1:4
%    subplot(2,2,k);
%    ylim([rhoMin, rhoMax]);
% end


% 
% 
% figure(4)
% for k = 1:4
%    subplot(2, 2, k)
%    %xlim([minCor(2), maxCor(2)])
%    ylim([minCor(1), maxCor(1)])
%    axis square
% end




% P = P1r;
cnt = 0;
rhoMin = NaN;
rhoMax = NaN;
for set =  ["GPS", "GAL", "GLO", "BEI"]
% rhoMin = min(rhoMin, min(min(P.RHO.(set))));
% rhoMax = max(rhoMax, max(max(P.RHO.(set))));
cnt = cnt + 1;
    
time = length(P.RHO.(set)(1,:));
xTicks = 0 : 900 : time;

sat = zeros(1,time);
for s = 1 : time
    sat(s) = sum(not(isnan(P.RHO.(set)(:,s))));
end

%% WLMS
L = length(P.RHO.(set)(:,1));
sigmaUERE = zeros(1, L);
for c = 1 : L
    r = P.RHO.(set)(c,:);
    r = r(~(isnan(r)));
    diffR = diff(r, 2);
    sigmaUERE(c)=std(diffR);
end
R = diag(sigmaUERE).^2;

for n = 1 : time  % time instants
    
    % Initialization.
    xHat = [0, 0, 0, 0]; 
    
    % Variables
    rho = P.RHO.(set)(:,n);
    visSat = find(not(isnan(P.RHO.(set)(:,n))));  % Visible sat. indexes
    H = zeros(length(visSat),4);
    H(:,4) = 1;
                         
    for k = 1:7  % Iterations of algorithm
        
        % Step 0: rhoHat and H_matrix evaluation
        for s = 1 : length(visSat)  % Iterations on satellites
            e = P.SAT_POS_ECEF.(set)(visSat(s)).pos(n,:) - xHat(1:3);
            rhoHat(s) = sqrt(sum(e.^2));
            H(s,1:3) = e/rhoHat(s);
        end  
        
        % Step 1: DeltaRhoHat  
        eRhoHat = rhoHat' - rho(visSat);
  
        % Step 2: H_w
        W = inv(R(visSat, visSat));
        H_w = inv(H.'*W*H)*(H.')*W;
        
        % Step 3: DeltaXHat
        eXHat = H_w * eRhoHat; 
       
        % Step 4: xHat
        xHat = xHat + eXHat.';
                
    end
    
    H_w = H_w.';
    G = inv(H_w.'*H_w);
    pdop(n) = sqrt(trace(G(1:end-1, 1:end-1)));
    gdop(n) = sqrt(trace(G));
    sigmaX(n) = sqrt(trace(inv(H.'*H) * H.'*R(visSat, visSat)*H*inv(H.'*H)));
    
    pos(n,1:4) = xHat;  %pos è in x,y,z
    clear rhoHat
    clear xHat
end

fprintf("WLMS sigmaX for %s =\n", [set])
% display(posError')
% display(sqrt(sum(posError.^2))')

for tPos = 1 : length(pos)
    coord(tPos, 1:3) = ecef2lla(pos(tPos, 1:3));  %coord è in lon,lat,alt
end

mPdop(cnt,2) = mean(pdop);
cov_pos = cov(pos);
cov_pos = cov_pos(1:3, 1:3);
cov_pos = trace(cov_pos);
wlmsSigmaX(cnt) = sqrt(cov_pos);

%% Plot: Visible satellites vs time.
% figure(11)
% plot(sat, 'linewidth', 2)
% xlabel("Time")
% ylabel("Number of visible satellites")
% title("Number of visible satellites vs Time")
% ylim([min(sat)-1, max(sat)+1])
% xlim([0,time])
% xticks(xTicks);
% grid on

%% Plot: Pseudoranges vs time.
figure(12)
subplot(2, 2, cnt)
plot(P.RHO.(set)', 'linewidth', 1)
% plot([1:1:time], P.RHO.(set)', '.')
xlabel("Time")
ylabel("Pseudorange [m]")
title([set])
xlim([0,time])
xticks(xTicks);
grid on
grid minor


%% Plot: World Map Plot
figure(13)
hold on
WorldMapPlotter(coord(time,1), coord(time,2), place, 'World')

%% Plot: Scatter plot of positions
% Real position is the mean of all positions.
% realPos = [mean(pos(:,1)), mean(pos(:,2)), mean(pos(:,3))];
% realPosLla = ecef2lla(realPos);
% errors = realPos - pos(:,1:3);
% sigmaX = sqrt(trace(cov(errors)));

fprintf("Avg. PDOP = %.6f metres\n", mean(pdop))
fprintf("Avg. GDOP = %.6f metres\n", mean(gdop))
% fprintf("Error: %.2d metres\n", sigmaX);
fprintf("Avg. sigmaUERE = %.2d\n", nanmean(sigmaUERE));
mean(pdop)*(nanmean(sigmaUERE));

figure(4)
subplot(2, 4, cnt+4)
hold on
geoshow(coord(:,1), coord(:,2),'DisplayType', 'point', 'Marker', '.',...
        'MarkerEdgeColor', [0, 0.4470, 0.7410])
geoshow(realPosLla(1), realPosLla(2),'DisplayType', 'point',...
        'MarkerEdgeColor', 'black', 'Marker','o',...
        'MarkerFaceColor', [0.8500, 0.3250, 0.0980])
legend('Estimated', 'Real')
xlabel('Longitude')
ylabel('Latitude')
title(set)%'Positions')
axis equal; axis square; grid on;


%% Plot: PDOP vs time.
figure(15)
subplot(2, 2, cnt)
hold on
plot(gdop)
plot(sigmaX)
title(set)%"Dilution of precision")
legend('GDOP', '\sigma_X', 'location', 'best')
xlabel('Time')
% ylabel('PDOP')
xlim([0,time])
xticks(xTicks);
grid on
grid minor



%% end
std(errors)
clear R  sigmaUERE sigmaX H G 
end


figure(4)
linkaxes([subplot(2,4,1),subplot(2,4,2),subplot(2,4,3),subplot(2,4,4),...
    subplot(2,4,5),subplot(2,4,6),subplot(2,4,7),subplot(2,4,8)],'xy')
%linkaxes([subplot(2,4,1),subplot(2,4,2),subplot(2,4,3),...
%    subplot(2,4,5),subplot(2,4,6),subplot(2,4,7)],'xy')

% figure(14)
% linkaxes([subplot(2,2,1),subplot(2,2,2),subplot(2,2,3),subplot(2,2,4)],'xy')

