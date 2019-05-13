close all;
clear;
clc;
format long;

a = [6376985;
     6377276;
     6377397;
     6378338;
     6378249;
     6378140;
     6378388;  % Hayford
     6378245;  % Krassovsky
     6378137];  % WGS84
 
 alpha = [1/308.6;
          1/300.8;
          1/299.1528128;
          1/288.5;
          1/293.5;
          1/298.3;
          1/297.0;
          1/298.3;
          1/298.257223563];
      
c = a .* (1 - alpha);
e2 = (a.^2 - c.^2) ./ (a.^2);
e_2 = (a.^2 - c.^2) ./ (c.^2);


%% First set of coordinates
% First task.
 k = 1;
 for i = [7 9]
    
    lat = deg2rad(dms2degrees([44 45 01.03930])); % Phi
    lon = deg2rad(dms2degrees([7 24 29.20335]));  % Lambda
    h = 322.4909;
        
    W1(k,1) = sqrt(1 - e2(i)*(sin(lat))^2);
    X1(k,1) = (a(i)/W1(k,1) + h) .* cos(lat) *cos(lon);
    Y1(k,1) = (a(i)/W1(k,1) + h) .* cos(lat) * sin(lon);
    Z1(k,1) = ((a(i)/W1(k,1)) * (1 - e2(i)) + h) * sin(lat);
    
    k = k + 1;
 end
  
 % Second task (+2000).
 k = 1;
 for i = [7 9]      
    
    lat = deg2rad(dms2degrees([44 45 01.03930]));  %Phi
    lon = deg2rad(dms2degrees([7 24 29.20335]));  %Lambda
    h = 322.4909;
    h = h + 2000;
        
    W1(k,2) = sqrt(1 - e2(i)*(sin(lat))^2);
    X1(k,2) = (a(i)/W1(k,2) + h) .* cos(lat) *cos(lon);
    Y1(k,2) = (a(i)/W1(k,2) + h) .* cos(lat) * sin(lon);
    Z1(k,2) = ((a(i)/W1(k,2)) * (1 - e2(i)) + h) * sin(lat);
    
    
    k = k + 1;
 end  

rowNames = {'Hayford','WGS84'};
colNames = {'Normal','m2000'};
X1 = array2table(X1,'RowNames',rowNames,'VariableNames',colNames)
Y1 = array2table(Y1,'RowNames',rowNames,'VariableNames',colNames)
Z1 = array2table(Z1,'RowNames',rowNames,'VariableNames',colNames)
 
 
%% Second set of coordinates
% First task.
k = 1;
  for i = [7 9]
      
    lat = deg2rad(dms2degrees([44 47 10.90505]));  % Phi
    lon = deg2rad(dms2degrees([7 30 26.53939]));  % Lambda
    h = 305.7367;
        
    W2(k,1) = sqrt(1-e2(i)*(sin(lat))^2);
    X2(k,1) = (a(i)/W2(k,1)+h).*cos(lat)*cos(lon);
    Y2(k,1) = (a(i)/W2(k,1)+h).*cos(lat)*sin(lon);
    Z2(k,1) = ((a(i)/W2(k,1))*(1-e2(i))+h)*sin(lat);
    
    k = k + 1;
  end
  
  % Second task (+2000).
k = 1;
  for i = [7 9]
    lat = dms2degrees([44 47 10.90505]); % Phi
    lon = dms2degrees([7 30 26.53939]);  % Lambda
    h = 2305.7367;
        
    W2(k,2) = sqrt(1-e2(i)*(sin(lat))^2);
    X2(k,2) = (a(i)/W2(k,2)+h).*cos(lat)*cos(lon);
    Y2(k,2) = (a(i)/W2(k,2)+h).*cos(lat)*sin(lon);
    Z2(k,2) = ((a(i)/W2(k,2))*(1-e2(i))+h)*sin(lat);
    
    k = k + 1;
  end
  
X2 = array2table(X2,'RowNames',rowNames,'VariableNames',colNames)
Y2 = array2table(Y2,'RowNames',rowNames,'VariableNames',colNames)
Z2 = array2table(Z2,'RowNames',rowNames,'VariableNames',colNames)
    
%hayford_1 = [X1(1) X11(1); Y1(1) Y11(1); Z1(1) Z11(1)];
%WGS84_1 = [X1(2) X11(2); Y1(2) Y11(2); Z1(2) Z11(2)];