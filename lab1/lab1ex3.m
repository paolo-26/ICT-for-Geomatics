clc; clear; close all;
format

filename2="dataset_HEL_20180329T160900.mat"

P2 =load (filename2)
%place='    TRN'
time = length(P2.RHO.GPS(1,:));
xTicks = [0:400:3600];
for s = 1:time
    sat(s) = sum(not(isnan(P2.RHO.GPS(:,s))));
end

figure(1)
plot(sat, 'linewidth', 2)
xlabel("Time")
ylabel("Number")
title("Number of visible satellites vs Time")
ylim([11.5,17.5])
grid on
figure(2)

plot([1:1:time],P2.RHO.GPS', '.')
xlabel("Time")
ylabel("Pseudorange")
title("Pseudoranges vs Time")
xticks(xTicks);
grid on

%initialize
%x_hat = zeros(1,4)
x_hat=[0,0,0,1]; 

for n=1:32
    vars(n) = nanvar(P2.RHO.GPS(n,:));
end

m=1;
satells = [5:14]
for maxSat = satells

    %indici= find(not(isnan(vars)));
    %indici(maxSat+1:end)=[];
    
    ml=0; %solo per il grafico

        for n = 1:time  %istanti di tempo
            x_hat=[0,0,0,1]; 
            rho = P2.RHO.GPS(:,n);

            n_sat = sat(n); %numero di satelli visibili al tempo n

            visSat = find(not(isnan(P2.RHO.GPS(:,n)))); %indici dei satelliti visibili
            if (length(visSat)>ml)
                ml=length(visSat);
            end
            
            if (length(visSat)>maxSat)
                visSat(maxSat+1:end)=[];
                n_sat = length(visSat);
            end
            %if (n_sat~=length(visSat))
            %    disp('Error');
            %end

            % H matrix initialization.
            H = zeros(length(visSat),4);
            H(:,4) = 1;     

            % R matrix
            R = zeros(length(visSat), length(visSat));
            for i = 1:length(visSat)
                R(i,i) = vars(visSat(i));
            end
            
            % W matrix.
            W = inv(R);
    
            for k = 1:12  % iterazioni all'instante n

                for s = 1:length(visSat) % iterazioni su ogni satellite

                    e=P2.SAT_POS_ECEF.GPS(visSat(s)).pos(n,:)-x_hat(1:3);
                    rho_hat(s)=sqrt(sum(e.^2));
                    H(s,1:3)=[e/rho_hat(s)];
                    H(:,4)=1; 

                end

                H_w = inv(H'*W*H)*H'*W;
                e_rho_hat = rho_hat' - rho(visSat);
                e_x_hat = H_w*e_rho_hat;
                x_hat=x_hat+e_x_hat';

            end
%             H2=H(1:end-1,1:end-1);
%             Q=inv(H2'*H2);
%             pdop(n) = sqrt(trace(Q));
            Q = inv(H'*H);
            pdop(n) = sqrt(trace(Q(1:end-1, 1:end-1)));

            pos(n,1:4) = x_hat;  %pos in (x, y, z)
            clear rho_hat
        end

    for n=1:length(pos)
        coord(n,1:3) = ecef2lla(pos(n,1:3));  %coord in (lon, lat, alt)
    end

    %figure

    realPos = [mean(pos(:,1)), mean(pos(:,2)), mean(pos(:,3))];
    realPosLla=ecef2lla(realPos);
    errors = realPos - pos(:,1:3);
    %plot(errors)

    %H2=H(1:end-1,1:end-1);
    %pdop = sqrt(trace(inv(H2'*H2)));

    errCovMat(m) = sqrt(trace(cov(errors)));

    meanPdop(m) = mean(pdop);
    
    % Print results.
    fprintf("No. of satellites: %d\n", maxSat)
    fprintf("  PDOP  = %.2d metres\n", meanPdop(m))
    fprintf("  Error = %.2d metres\n\n", errCovMat(m))
    
    % Save axes limits for plots.
    if (m == 1)
       lim1 = [min(coord(:,1)), max(coord(:,1))];
       lim2 = [min(coord(:,2)), max(coord(:,2))];
    end
        
    m = m + 1;

figure
%subplot(1,2,1)
%WorldMapPlotter(realPosLla(1),realPosLla(2), place, 'Italy')
%subplot(1,2,2)
hold on
geoshow(coord(:,1), coord(:,2),'DisplayType', 'point', 'Marker', '.')
geoshow(realPosLla(1), realPosLla(2),'DisplayType', 'point',...
                                            'MarkerEdgeColor','blue',...
                                            'Marker', 'o',...
                                            'MarkerFaceColor', 'green')
legend('Estimated positions', 'Real position')
xlabel('Longitude')
ylabel('Latitude')
title([num2str(maxSat) ' satelliti'])

ylim(lim1)
xlim(lim2)
grid on

end

figure
hold on

k = 1;
for n = 1 : 32
    d(n,:) = diff(diff(P2.RHO.GPS(n,:)));
    
    if (not(isnan(d(n,1))))
        subplot(ml,1,k)
        plot(d(n,:))
        k=k+1;
    end
end

% figure(10)
% subplot(1,2,1)
% plot(diff(diff(P2.RHO.GPS(2,:))))
% 
% subplot(1,2,2)
% plot(diff(P2.RHO.GPS(2,:)))

figure
plot(pdop)
title("Position dilution of precision")
grid minor
xlabel('Time')
ylabel('pdop')
xlim([0,3600])
xticks(xTicks);

figure
hold on
plot(satells, meanPdop)
plot(satells, errCovMat)
legend('PDOP', 'COV ERR')
xlabel("Number of satellites")
grid on

