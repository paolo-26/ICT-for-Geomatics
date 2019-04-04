% Code for tasks 1-2

clc; clear; close all;
format short

% NominalUEREfixed.
filename = "dataset_HEL_20180328T122158";  % Helsinki
place = '     HEL';  % Label for world plot

P1 = load(filename);
time = length(P1.RHO.GPS(1,:));
xTicks = [0 : 400 : time];

sat = zeros(1,time);
for s = 1 : time
    sat(s) = sum(not(isnan(P1.RHO.GPS(:,s))));
end

%% Recursive Least Mean Square

for n = 1 : time  % time instants
    
    % Initialization.
    xHat = [0, 0, 0, 1]; 
    
    % Variables
    rho = P1.RHO.GPS(:,n);
    visSat = find(not(isnan(P1.RHO.GPS(:,n))));  %visible satellites indexes
    H = zeros(length(visSat),4);
    H(:,4) = 1;
    
    % Remove some satellites.
    %visSat(7:end)=[];
            %satNo = length(visSat);

            %if (satNo~=length(visSat))
            %    disp('Error');
            %end
       
                     
    for k = 1:10  % Iterations of algorithm
        
        % Step 0: rhoHat and H_matrix evaluation
        for s = 1 : length(visSat)  % Iterations on satellites
            e = P1.SAT_POS_ECEF.GPS(visSat(s)).pos(n,:) - xHat(1:3);
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
    
    Q = inv(H'*H);
    pdop(n) = sqrt(trace(Q(1:end-1, 1:end-1)));
    
    pos(n,1:4) = xHat;  %pos è in x,y,z
    clear rhoHat
    clear xHat
end

for tPos = 1 : length(pos)
    coord(tPos, 1:3) = ecef2lla(pos(tPos, 1:3));  %coord è in lon,lat,alt
end

r = P1.RHO.GPS(1,:);
r = r(~(isnan(r(1,:))));
diffR = diff(r, 2);
mean(diffR)

r = P1.RHO.GPS(12,:);
r = r(~(isnan(r(1,:))));
diffR = diff(r, 2);
mean(diffR)

% Save position for Google Earth
filename = 'Helsinki.m';
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
ylim([11.5, 17.5])
xlim([0,time])
xticks(xTicks);
grid on

%% Plot: Pseudoranges vs time.
figure(2)
plot([1:1:time], P1.RHO.GPS', '.')
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
err = sqrt(trace(cov(errors)));
fprintf("PDOP = %.2d metres\n", pdop)
fprintf("Error: %.2d metres\n", err)

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
