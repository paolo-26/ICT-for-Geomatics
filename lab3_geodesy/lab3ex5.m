close all;
clear;
clc;
format long;

datumA = [...
4472675.1340 677732.7296 4481406.6192;  %ALES      
4478841.8179 645666.0585 4480078.0380;  %ASTI      
4429584.8191 626326.2574 4531541.7418;  %BIEL      
4446830.5758 601969.5185 4517928.9157;  %CAST      
4525690.8631 600123.4090 4439976.2880;  %CUNE      
4474374.0522 550698.8143 4497968.8137;  %GRAV      
4523227.3918 621792.1742 4439482.9654;  %MOND      
4431899.3484 671366.9963 4522512.0588;  %NOVA      
4444603.4927 714785.8615 4503373.0213;  %PAVI      
4504896.3605 605936.8257 4459842.1457;  %SAVI      
4472544.5785 601634.0764 4492545.0107;  %TORI      
4443192.0735 657757.3082 4513434.4540];  %VERC  
datumB = [...
4472674.8193 677733.0708 4481406.8974;  %ALES      
4478841.5110 645666.3996 4480078.3081;  %ASTI      
4429584.5150 626326.5832 4531541.9794;  %BIEL      
4446830.2784 601969.8594 4517929.1541;  %CAST      
4525690.5486 600123.7968 4439976.5664;  %CUNE      
4474373.7679 550699.1715 4497969.0384;  %GRAV      
4523227.0920 621792.5562 4439483.2501;  %MOND      
4431899.0641 671367.3148 4522512.3138;  %NOVA      
4444603.1685 714786.1790 4503373.2881;  %PAVI      
4504896.0640 605937.2100 4459842.4154;  %SAVI      
4472544.2920 601634.4325 4492545.2622;  %TORI      
4443191.7627 657757.6360 4513434.7201];  %VERC  

for i = 1 : size(datumA, 1)
    Xa = datumA(i, 1) / 1e6;
    Ya = datumA(i, 2) / 1e6;
    Za = datumA(i, 3) / 1e6;
    Xb = datumB(i, 1) / 1e6;
    Yb = datumB(i, 2) / 1e6;
    Zb = datumB(i, 3) / 1e6;
    
    newA = [1 0 0  Xa   0   -Za     Ya;
            0 1 0  Ya  -Za   0     -Xa;  
            0 0 1  Za   Ya   Xa     0 ];
        
    loNew = [Xb-Xa; Yb-Ya; Zb-Za]; 
    
    if i == 1
        A = newA;
        lo = loNew;
    else
        A = [A; newA];
        lo = [lo; loNew];
    end
    
end

x = (A.'*A)\(A.'*lo);
x(1:3, 1) = x(1:3,1) * 1e6;
v = (A*x - lo) ;

Names = {'Tx (translation)', 'Ty (translation)', 'Tz (translation)', 'lambda (scaling)', 'Rx (rotation)', 'Ry (rotation)', 'Rz (rotation)'};
x = array2table(x, 'rowNames', Names)