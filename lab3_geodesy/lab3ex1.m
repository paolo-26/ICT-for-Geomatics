close all;
clear;
clc;
format shorte;

a = [6376985; 6377276; 6377397; 6378338; 6378249;
     6378140; 6378388; 6378245; 6378137];
 
alpha = [1/308.6; 1/300.8; 1/299.1528128; 1/288.5;
         1/293.5; 1/298.3; 1/297.0; 1/298.3;
         1/298.257223563];
      
c = a .* (1 - alpha);
e2 = (a.^2 - c.^2) ./ (a.^2);
e_2 = (a.^2 - c.^2) ./ (c.^2);

colNames = {'a', 'alpha', 'c', 'e2', 'e_2'};
rowNames = {'Delambre', 'Everest', 'Bessel', 'Fisher', 'Clarke',...
            'Helmert', 'Hayford', 'Krassovsky', 'WGS84'};
array2table([a alpha c e2 e_2], 'RowNames',rowNames,...
            'VariableNames', colNames)
