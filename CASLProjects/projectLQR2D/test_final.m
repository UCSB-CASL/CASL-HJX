% LQR2D_Visualization.m
% Creates a 6x3 subplot grid visualizing LQR value functions and errors
% with correct folder structure and domain range

clear all; close all; clc;

% === Configuration ===
gridSize = 80;           % Grid size to visualize
tFinal = 5;              % Final time
times = 0:5;             % Time points to visualize (t=0,1,2,3,4,5)
domainRange = 2;         % Domain range [-domainRange, domainRange]²
dt = 0.01;               % Time step dt value to visualize

% Output directory with proper grid size subfolder
output_dir = sprintf('./LQR2D_Output_dt%.2f/LQR2D_%d/', dt, gridSize);

% === Setup figure with LaTeX formatting ===
figure('Position', [50, 50, 420, 1200]);
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');

% === Create coordinate grid ===
[X1, X2] = meshgrid(linspace(-domainRange, domainRange, gridSize));

% === Loop over time points (t=0,1,2,3,4,5) ===
for timeIdx = 1:length(times)
    time = times(timeIdx);
    tau = tFinal - time;  % Convert to tau (goes backward)
    
    % === Load numerical solution ===
    phi_file_name = sprintf('phi_t%d.dat', time);
    phi_dir = fullfile(output_dir, 'phi', phi_file_name);
    
    try
        phi_file = readmatrix(phi_dir);
    catch
        warning('File not found: %s. Using zeros instead.', phi_dir);
        phi_file = zeros(gridSize);
    end
    
    % === Calculate exact solution ===
    V_exact = exactSolution(X1, X2, time, tFinal);
    
    % === Calculate error ===
    error = abs(V_exact - phi_file);
    max_error = max(error(:));
    l2_error = sqrt(mean(error(:).^2));
    
    % === Row index for this time step ===
    row = time + 1;
    
    % === Column 1: Exact Solution ===
    ax1 = subplot(6, 2, 2*(row-1) + 1);
    surf(X1, X2, V_exact, 'EdgeColor', 'none');
    title(sprintf('$\\textrm{Exact},\\ \\tau = %d$', tau), 'FontWeight', 'bold');
    xlabel('$x_1$', 'FontWeight', 'bold');
    ylabel('$x_2$', 'FontWeight', 'bold');
    zlabel('$V(x,\\tau)$', 'FontWeight', 'bold');
    view([-30, 30]);
    colormap(ax1, parula);
    grid off;
    set(ax1, 'LineWidth', 1.5);
    
    % === Column 2: Numerical Solution ===
    ax2 = subplot(6, 2, 2*(row-1) + 2);
    surf(X1, X2, phi_file, 'EdgeColor', 'none');
    title(sprintf('$\\textrm{Numerical},\\ \\tau = %d$', tau), 'FontWeight', 'bold');
    xlabel('$x_1$', 'FontWeight', 'bold');
    ylabel('$x_2$', 'FontWeight', 'bold');
    zlabel('$V(x,\\tau)$', 'FontWeight', 'bold');
    view([-30, 30]);
    colormap(ax2, parula);
    grid off;
    set(ax2, 'LineWidth', 1.5);
    
    % % === Column 3: Error Visualization ===
    % ax3 = subplot(6, 2, 2*(row-1) + 3);
    % surf(X1, X2, error, 'EdgeColor', 'none');
    % title(sprintf('$\\textrm{Error},\\ \\max = %.2e,\\ L_2 = %.2e$', max_error, l2_error), 'FontWeight', 'bold');
    % xlabel('$x_1$', 'FontWeight', 'bold');
    % ylabel('$x_2$', 'FontWeight', 'bold');
    % zlabel('$|\\textrm{Error}|$', 'FontWeight', 'bold');
    % view([-30, 30]);
    % colormap(ax3, hot);
    % grid off;
    % set(ax3, 'LineWidth', 1.5);
end

% === Add overall title with LaTeX ===
% sgtitle(sprintf('\\textbf{LQR2D Solution: Grid Size = %d, dt = %.3f, Domain [-2,2]$^2$}', gridSize, dt), 'FontSize', 16, 'Interpreter', 'latex');

% === Save figure in high quality ===
set(gcf, 'Color', 'white');
exportgraphics(gcf, sprintf('LQR2D_Visualization_Grid%d_dt%.3f.png', gridSize, dt), 'Resolution', 300);

% === Function for exact solution ===
function V_exact = exactSolution(X1, X2, t, tFinal)
    % Riccati-based exact solution from backward integration
    persistent K_sol t_sol
    if isempty(K_sol)
        A = [0 1; 0 0];
        B = [0; 1];
        Q = eye(2);
        R = eye(1);
        P_T = eye(2);  % terminal cost
        K0 = reshape(P_T, [], 1);
        opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
        [t_sol, K_sol] = ode45(@(t, K) ode_gain(t, K, A, B, Q, R), ...
                               [tFinal:-0.01:0], K0, opts);
    end
    tau = tFinal - t;
    K_interp = reshape(interp1(t_sol, K_sol, tau, 'pchip'), 2, 2);
    V_exact = 0.5 * (K_interp(1,1)*X1.^2 + 2*K_interp(1,2)*X1.*X2 + K_interp(2,2)*X2.^2);
end

function dK = ode_gain(~, K, A, B, Q, R)
    K_mat = reshape(K, size(A));
    dK_mat = -(A' * K_mat + K_mat * A - K_mat * B * (R \ B') * K_mat + Q);
    dK = reshape(dK_mat, [], 1);
end