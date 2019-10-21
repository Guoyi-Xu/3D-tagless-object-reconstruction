close all
clear
clc

% Draw the grid points.
figure;
GridPos = [0.6, 0.6; 0.6, 1.2; 0.6, 1.8; 0.6, 2.4; 0.6, 3.0;
    1.2, 0.6; 1.2, 1.2; 1.2, 1.8; 1.2, 2.4; 1.2, 3.0;
    1.8, 0.6; 1.8, 1.2; 1.8, 1.8; 1.8, 2.4; 1.8, 3.0;
    2.4, 0.6; 2.4, 1.2; 2.4, 1.8; 2.4, 2.4; 2.4, 3.0;
    3.0, 0.6; 3.0, 1.2; 3.0, 1.8; 3.0, 2.4; 3.0, 3.0];
scatter(GridPos(:, 1), GridPos(:, 2), 40, 'r', 'filled');
hold on

ResPos = [0.6, 0.72; 0.48, 0.84; 0.84, 1.86; 3.18, 2.76; 1.32, 3.4;
    0.9, 0.84; 1.35, 1.26; 1.19, 1.85; 0.72, 3.2; 1.18, 2.28;
    1.92, 1.05; 1.8, 1.04; 1.8, 1.8; 1.86, 2.28; 1.89, 2.28;
    2.25, 0.9; 2.31, 1.2; 2.34, 1.78; 2.34, 2.76; 2.43, 2.79;
    3.6, 0.42; 3.09, 1.26; 2.91, 2.34; 2.88, 2.48; 3.44, 2.08];
scatter(ResPos(:, 1), ResPos(:, 2), 40, 'b', 'filled');

legend('Grid Positions', 'Result Positions');

xlim([0, 3.8]);
ylim([0, 3.8]);
xlabel('x position (m)');
ylabel('y position (m)');
set(gca, 'fontweight', 'bold');

CalDist = sqrt((GridPos(:, 1) - ResPos(:, 1)).^2 + (GridPos(:, 2) - ResPos(:, 2)).^2);
CalDist = (reshape(CalDist, [5, 5]))';

figure;
GridPosGood = [0.6, 0.6; 0.6, 1.2; 0.6, 1.8;
    1.2, 0.6; 1.2, 1.2; 1.2, 1.8;
    1.8, 0.6; 1.8, 1.2; 1.8, 1.8; 1.8, 2.4;
    2.4, 0.6; 2.4, 1.2; 2.4, 1.8; 2.4, 2.4; 2.4, 3.0;
    3.0, 1.2; 3.0, 1.8; 3.0, 2.4];
scatter(GridPosGood(:, 1), GridPosGood(:, 2), 40, 'm', 'filled');
hold on

GridPosNotGood = [0.6, 2.4; 0.6, 3.0;
    1.2, 2.4; 1.2, 3.0;
    1.8, 3.0;
    
    3.0, 0.6; 3.0, 3.0];
scatter(GridPosNotGood(:, 1), GridPosNotGood(:, 2), 40, 'g', 'filled');

legend('Accuracy < 60 cm', 'Poor Accuracy');

xlim([0, 3.8]);
ylim([0, 3.8]);
xlabel('x position (m)');
ylabel('y position (m)');
set(gca, 'fontweight', 'bold');