close all; clear all; clc;

phi = deg2rad([0 30 37 45 60]);
lambda = deg2rad([0 : 0.5 : 3]);

figure(1)
hold on
grid on
for k = 1 : length(phi)
    m = 0.9996*(1 + (((lambda.^2)./2).*(cos(phi(k))).^2));

    plot([0 : 0.5 : 3], m, '.-', 'linewidth',2, 'Markersize', 20)
    if k == 1
        yticks(m)
    end
end

legend("equator", "case 1", "case 2", "case 3", "case 4",...
       'location', 'best');
ylabel("\it m_l");
xlabel("longitude");
xtickformat('degrees');
grid minor
title('Modulus of linear deformation vs longitude')