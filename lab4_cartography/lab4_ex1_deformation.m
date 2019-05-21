close all; clear all; clc;

phi = deg2rad([30 37 45 60]);
lambda = deg2rad([0 : 0.5 : 3]);

figure(1)
hold on
grid on
for k = 1 : length(phi)
    m = 0.9996*(1 + (((lambda.^2)./2).*(cos(phi(k))).^2));
    plot([0 : 0.5 : 3], m, 'o-', 'linewidth',2)
end

legend("case 1", "case 2", "case 3", "case 4", 'location','northwest');
ylabel("\it m_l");
yticks([0.9996 : 0.0001 : 1.0008])
xlabel("longitude");
xtickformat('degrees');
grid minor
title('Modulus of linear deformation vs longitude')