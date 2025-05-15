% LQR2D_Convergence_Analysis_Simple.m
% Simple convergence analysis for LQR2D problem using the actual results

clear all; close all; clc;

% === Configuration ===
gridSizes = [20, 40, 80, 160, 320];  % Grid sizes to analyze
dtValues = [0.1, 0.01, 0.001];       % Time step values
tFinal = 5;                          % Final time
times = 0:5;                         % Time points (t=0,1,2,3,4,5)
tauValues = tFinal - times;          % τ values (τ=5,4,3,2,1,0)

% === Hard-coded actual results (from your output) ===
% Format: maxErrors(grid_idx, time_idx, dt_idx)
maxErrors = zeros(length(gridSizes), length(times), length(dtValues));
l2Errors = zeros(length(gridSizes), length(times), length(dtValues));

% dt = 0.1 results
maxErrors(1, 1, 1) = 4.792244e-06; l2Errors(1, 1, 1) = 2.174573e-06;
maxErrors(1, 2, 1) = 5.074858e-01; l2Errors(1, 2, 1) = 2.362350e-01;
maxErrors(1, 3, 1) = 1.067711e-01; l2Errors(1, 3, 1) = 4.967369e-02;
maxErrors(1, 4, 1) = 2.914019e-03; l2Errors(1, 4, 1) = 1.348099e-03;
maxErrors(1, 5, 1) = 3.236888e+01; l2Errors(1, 5, 1) = 1.635179e+00;
maxErrors(1, 6, 1) = NaN;          l2Errors(1, 6, 1) = NaN;

maxErrors(2, 1, 1) = 4.996713e-06; l2Errors(2, 1, 1) = 2.414684e-06;
maxErrors(2, 2, 1) = 5.102691e-01; l2Errors(2, 2, 1) = 2.252799e-01;
maxErrors(2, 3:6, 1) = NaN;        l2Errors(2, 3:6, 1) = NaN;

maxErrors(3, 1, 1) = 4.992790e-06; l2Errors(3, 1, 1) = 2.283187e-06;
maxErrors(3, 2, 1) = NaN;          l2Errors(3, 2, 1) = NaN;
maxErrors(3, 3:6, 1) = NaN;        l2Errors(3, 3:6, 1) = NaN;

maxErrors(4, 1, 1) = 4.999802e-06; l2Errors(4, 1, 1) = 2.281913e-06;
maxErrors(4, 2, 1) = NaN;          l2Errors(4, 2, 1) = NaN;
maxErrors(4, 3:6, 1) = NaN;        l2Errors(4, 3:6, 1) = NaN;

maxErrors(5, 1, 1) = 4.999853e-06; l2Errors(5, 1, 1) = 2.249841e-06;
maxErrors(5, 2, 1) = NaN;          l2Errors(5, 2, 1) = NaN;
maxErrors(5, 3:6, 1) = NaN;        l2Errors(5, 3:6, 1) = NaN;

% dt = 0.01 results
maxErrors(1, 1, 2) = 4.792244e-06; l2Errors(1, 1, 2) = 2.174573e-06;
maxErrors(1, 2, 2) = 5.186658e-01; l2Errors(1, 2, 2) = 2.373839e-01;
maxErrors(1, 3, 2) = 1.080676e-01; l2Errors(1, 3, 2) = 5.021220e-02;
maxErrors(1, 4, 2) = 2.904450e-03; l2Errors(1, 4, 2) = 1.353927e-03;
maxErrors(1, 5, 2) = 2.750804e-03; l2Errors(1, 5, 2) = 1.279235e-03;
maxErrors(1, 6, 2) = 6.195513e-04; l2Errors(1, 6, 2) = 2.864151e-04;

maxErrors(2, 1, 2) = 4.996713e-06; l2Errors(2, 1, 2) = 2.414684e-06;
maxErrors(2, 2, 2) = 5.191935e-01; l2Errors(2, 2, 2) = 2.264302e-01;
maxErrors(2, 3, 2) = 1.082867e-01; l2Errors(2, 3, 2) = 4.789522e-02;
maxErrors(2, 4, 2) = 2.918363e-03; l2Errors(2, 4, 2) = 1.291606e-03;
maxErrors(2, 5, 2) = 2.753889e-03; l2Errors(2, 5, 2) = 1.219973e-03;
maxErrors(2, 6, 2) = 6.145625e-04; l2Errors(2, 6, 2) = 2.734281e-04;

maxErrors(3, 1, 2) = 4.992790e-06; l2Errors(3, 1, 2) = 2.283187e-06;
maxErrors(3, 2, 2) = 5.191819e-01; l2Errors(3, 2, 2) = 2.209928e-01;
maxErrors(3, 3, 2) = 1.083271e-01; l2Errors(3, 3, 2) = 4.674533e-02;
maxErrors(3, 4, 2) = 2.920240e-03; l2Errors(3, 4, 2) = 1.260437e-03;
maxErrors(3, 5, 2) = 2.758121e-03; l2Errors(3, 5, 2) = 1.190590e-03;
maxErrors(3, 6, 2) = 6.216523e-04; l2Errors(3, 6, 2) = 2.666771e-04;

maxErrors(4:5, :, 2) = NaN; l2Errors(4:5, :, 2) = NaN;

% dt = 0.001 results
maxErrors(1, 1, 3) = 4.792244e-06; l2Errors(1, 1, 3) = 2.174573e-06;
maxErrors(1, 2, 3) = 5.186658e-01; l2Errors(1, 2, 3) = 2.373839e-01;
maxErrors(1, 3, 3) = 1.079276e-01; l2Errors(1, 3, 3) = 5.015703e-02;
maxErrors(1, 4, 3) = 3.314019e-03; l2Errors(1, 4, 3) = 1.442257e-03;
maxErrors(1, 5, 3) = 2.822692e-03; l2Errors(1, 5, 3) = 1.285670e-03;
maxErrors(1, 6, 3) = 6.258556e-04; l2Errors(1, 6, 3) = 2.894285e-04;

maxErrors(2, 1, 3) = 4.996713e-06; l2Errors(2, 1, 3) = 2.414684e-06;
maxErrors(2, 2, 3) = 5.191935e-01; l2Errors(2, 2, 3) = 2.264301e-01;
maxErrors(2, 3, 3) = 1.081467e-01; l2Errors(2, 3, 3) = 4.784266e-02;
maxErrors(2, 4, 3) = 3.323123e-03; l2Errors(2, 4, 3) = 1.375644e-03;
maxErrors(2, 5, 3) = 2.818931e-03; l2Errors(2, 5, 3) = 1.226203e-03;
maxErrors(2, 6, 3) = 6.245625e-04; l2Errors(2, 6, 3) = 2.760971e-04;

maxErrors(3, 1, 3) = 4.992790e-06; l2Errors(3, 1, 3) = 2.283187e-06;
maxErrors(3, 2, 3) = 5.191819e-01; l2Errors(3, 2, 3) = 2.209927e-01;
maxErrors(3, 3, 3) = 1.081871e-01; l2Errors(3, 3, 3) = 4.669390e-02;
maxErrors(3, 4, 3) = 3.319592e-03; l2Errors(3, 4, 3) = 1.342411e-03;
maxErrors(3, 5, 3) = 2.824931e-03; l2Errors(3, 5, 3) = 1.196677e-03;
maxErrors(3, 6, 3) = 6.282320e-04; l2Errors(3, 6, 3) = 2.695377e-04;

maxErrors(4, 1, 3) = 4.999802e-06; l2Errors(4, 1, 3) = 2.281913e-06;
maxErrors(4, 2, 3) = 5.192095e-01; l2Errors(4, 2, 3) = 2.182857e-01;
maxErrors(4, 3, 3) = 1.082010e-01; l2Errors(4, 3, 3) = 4.612196e-02;
maxErrors(4, 4, 3) = 3.323033e-03; l2Errors(4, 4, 3) = 1.325927e-03;
maxErrors(4, 5, 3) = 2.823868e-03; l2Errors(4, 5, 3) = 1.182134e-03;
maxErrors(4, 6, 3) = 6.269436e-04; l2Errors(4, 6, 3) = 2.662222e-04;

maxErrors(5, 1, 3) = 4.999853e-06; l2Errors(5, 1, 3) = 2.249841e-06;
maxErrors(5, 2, 3) = 5.192019e-01; l2Errors(5, 2, 3) = 2.169351e-01;
maxErrors(5, 3, 3) = 1.082014e-01; l2Errors(5, 3, 3) = 4.583656e-02;
maxErrors(5, 4, 3) = 3.322708e-03; l2Errors(5, 4, 3) = 1.317772e-03;
maxErrors(5, 5, 3) = 2.822770e-03; l2Errors(5, 5, 3) = 1.174810e-03;
maxErrors(5, 6, 3) = 6.306307e-04; l2Errors(5, 6, 3) = 2.645248e-04;

% === Create a stability table ===
fprintf('=== Stability Analysis ===\n');
fprintf('Grid Size | dt=0.1 | dt=0.01 | dt=0.001\n');
fprintf('---------+--------+---------+--------\n');

for grid_idx = 1:length(gridSizes)
    gridSize = gridSizes(grid_idx);
    
    % Count stable time steps for each dt
    stableCount1 = sum(~isnan(maxErrors(grid_idx, :, 1)));
    stableCount2 = sum(~isnan(maxErrors(grid_idx, :, 2)));
    stableCount3 = sum(~isnan(maxErrors(grid_idx, :, 3)));
    
    fprintf('%7d  | %6d | %7d | %6d\n', gridSize, stableCount1, stableCount2, stableCount3);
end
fprintf('\n');

% === Create convergence plots ===
% Focus on dt=0.001 where all grid sizes work
dt_idx = 3;  % dt=0.001

% === Plot final time convergence ===
figure;
% Max error convergence
maxErrorsFinal = squeeze(maxErrors(:, 6, dt_idx));  % Final time (tau=0)
valid_idx = ~isnan(maxErrorsFinal);
loglog(gridSizes(valid_idx), maxErrorsFinal(valid_idx), 'o-', 'LineWidth', 2);
grid on;
title(sprintf('Maximum Error at Final Time (\\tau=0, dt=%.3f)', dtValues(dt_idx)));
xlabel('Grid Size (N)');
ylabel('Maximum Error');

% Add reference line for 2nd order convergence
hold on;
if sum(valid_idx) >= 2
    x_ref = [gridSizes(find(valid_idx, 1, 'first')), gridSizes(find(valid_idx, 1, 'last'))];
    y_ref = maxErrorsFinal(valid_idx(1)) * (x_ref(1)./x_ref).^2;
    loglog(x_ref, y_ref, 'k--', 'DisplayName', 'O(h^2) Reference');
    legend('Numerical Results', 'O(h^2) Reference');
end

print -dpng LQR2D_MaxError_Convergence.png

% Plot L2 error convergence
figure;
l2ErrorsFinal = squeeze(l2Errors(:, 6, dt_idx));  % Final time (tau=0)
valid_idx = ~isnan(l2ErrorsFinal);
loglog(gridSizes(valid_idx), l2ErrorsFinal(valid_idx), 'o-', 'LineWidth', 2);
grid on;
title(sprintf('L2 Error at Final Time (\\tau=0, dt=%.3f)', dtValues(dt_idx)));
xlabel('Grid Size (N)');
ylabel('L2 Error');

% Add reference line for 2nd order convergence
hold on;
if sum(valid_idx) >= 2
    x_ref = [gridSizes(find(valid_idx, 1, 'first')), gridSizes(find(valid_idx, 1, 'last'))];
    y_ref = l2ErrorsFinal(valid_idx(1)) * (x_ref(1)./x_ref).^2;
    loglog(x_ref, y_ref, 'k--', 'DisplayName', 'O(h^2) Reference');
    legend('Numerical Results', 'O(h^2) Reference');
end

print -dpng LQR2D_L2Error_Convergence.png

% === Create error evolution plot ===
figure;
% Pick a specific grid size
grid_idx = 3;  % Grid size 80
gridSize = gridSizes(grid_idx);

% Plot error evolution over time
semilogy(times, squeeze(maxErrors(grid_idx, :, dt_idx)), 'o-', 'LineWidth', 2, 'DisplayName', 'Max Error');
hold on;
semilogy(times, squeeze(l2Errors(grid_idx, :, dt_idx)), 's-', 'LineWidth', 2, 'DisplayName', 'L2 Error');
grid on;
% title(sprintf('Error Evolution (Grid=%d, dt=%.3f)', gridSize, dtValues(dt_idx)));
xlabel('Time (t)');
ylabel('Error');
legend('Location', 'best');

print -dpng LQR2D_Error_Evolution.png

% === Compare different dt values at final time ===
figure;
% Only use grid sizes up to 80 where we have data for dt=0.01
valid_grids = 1:3;  % Grid sizes 20, 40, 80

% Plot max error for each dt
for dt_idx = 2:3  % dt=0.01 and dt=0.001
    loglog(gridSizes(valid_grids), maxErrors(valid_grids, 6, dt_idx), 'o-', 'LineWidth', 2, ...
        'DisplayName', sprintf('dt=%.3f', dtValues(dt_idx)));
    hold on;
end
grid on;
title('Time Step Comparison: Maximum Error at Final Time');
xlabel('Grid Size (N)');
ylabel('Maximum Error');
legend('Location', 'best');

print -dpng LQR2D_TimeStep_Comparison.png

fprintf('Convergence analysis complete. Figures saved.\n');