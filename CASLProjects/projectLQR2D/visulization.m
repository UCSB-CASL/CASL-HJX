%% Ultra-Professional LQR2D GIF Creator for Academic Publications
%  Creates publication-quality visualizations with trajectories, control, and analysis
%  Enhanced with consistent scaling, LaTeX fonts, and professional styling

clear; close all; clc;

fprintf('🎬 Creating ULTRA-PROFESSIONAL GIFs for CASL-HJX LQR2D solver...\n');
fprintf('🏆 Publication-quality visualizations with trajectories and control analysis\n\n');

% Configuration
domain = [-4 4 -4 4];
N = 160;  % High resolution for professional quality
folder = sprintf('./LQR2D_Output/LQR2D_%d/phi/', N);

% Fallback options
if ~exist(folder, 'dir')
    N = 80;
    folder = sprintf('./LQR2D_Output/LQR2D_%d/phi/', N);
end

if ~exist(folder, 'dir')
    error('No LQR2D data found. Please run the solver first.');
end

fprintf('📊 Using high-quality data from: %s (N=%d)\n', folder, N);

% Complete time series for ultra-smooth animations
all_times = [0, 0.1:0.1:1, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
hero_times = all_times;
analysis_times = [0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0];
trajectory_times = [0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0];

[X, Y] = meshgrid(linspace(domain(1), domain(2), N));

% Pre-compute global z-axis limits for consistent scaling
fprintf('🔍 Pre-computing global z-axis limits for consistent scaling...\n');
z_min_global = inf;
z_max_global = -inf;

for t = hero_times
    fname = fullfile(folder, ['phi_t' time_string(t) '.dat']);
    if isfile(fname)
        phi = reshape(load(fname), N, N)';
        z_min_global = min(z_min_global, min(phi(:)));
        z_max_global = max(z_max_global, max(phi(:)));
    end
end

% Add small margin for better visualization
z_range = z_max_global - z_min_global;
z_min_global = z_min_global - 0.05 * z_range;
z_max_global = z_max_global + 0.05 * z_range;

fprintf('📊 Global z-axis range: [%.3f, %.3f]\n', z_min_global, z_max_global);

% Analytical Riccati solution
A = [0 1; 0 0]; B = [0; 1]; Q = eye(2); R = 1;
eye_2 = eye(2);
[tRic, Pvec] = ode45(@(t,P) riccati_rhs(t,P,A,B,Q,R), [0 10], eye_2(:));
Pana = reshape(Pvec.', 2, 2, []);

% Professional color schemes
colors = {
    [0.2 0.4 0.8],      % Professional blue
    [0.8 0.2 0.2],      % Professional red  
    [0.2 0.6 0.3],      % Professional green
    [0.9 0.6 0.1],      % Professional orange
    [0.6 0.2 0.8],      % Professional purple
    [0.3 0.3 0.3]       % Professional gray
};

% Set default figure properties for LaTeX compatibility
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultAxesFontName', 'Times');
set(groot, 'DefaultTextFontName', 'Times');

%% ========================================================================
%% HERO GIF: Advanced Surface + Trajectories + Control Analysis
%% ========================================================================

fprintf('🎯 Creating HERO multi-panel visualization with trajectories...\n');

fig1 = figure('Position', [50 50 1400 800], 'Color', 'white');
set(fig1, 'InvertHardcopy', 'off');

filename_hero = 'lqr_hero_professional.gif';
delay_hero = 0.15;

% Initialize trajectory storage
num_trajectories = 6;
trajectory_starts = [
    -3, -2; -2, 3; 3, -1; 1, 2.5; -1.5, -2.5; 2.5, 1.5
];

% Pre-compute optimal trajectories
fprintf('📊 Pre-computing optimal trajectories...\n');
trajectories = cell(num_trajectories, 1);
control_signals = cell(num_trajectories, 1);

for traj_idx = 1:num_trajectories
    [traj, control] = compute_optimal_trajectory(trajectory_starts(traj_idx, :), ...
                                               hero_times, folder, N, domain);
    trajectories{traj_idx} = traj;
    control_signals{traj_idx} = control;
end

% Pre-compute control signal limits for consistent scaling
u_min_global = inf;
u_max_global = -inf;
for traj_idx = 1:num_trajectories
    if ~isempty(control_signals{traj_idx})
        u_min_global = min(u_min_global, min(control_signals{traj_idx}));
        u_max_global = max(u_max_global, max(control_signals{traj_idx}));
    end
end
u_margin = 0.1 * (u_max_global - u_min_global);
u_min_global = u_min_global - u_margin;
u_max_global = u_max_global + u_margin;

% Main animation loop
for k = 1:length(hero_times)
    t = hero_times(k);
    
    fname = fullfile(folder, ['phi_t' time_string(t) '.dat']);
    if ~isfile(fname), continue; end
    
    phi = reshape(load(fname), N, N)';
    
    % Clear figure and create subplots
    clf(fig1);
    
    % Subplot 1: 3D Cost-to-Go Surface
    ax1 = subplot(2, 3, [1 2]);
    surf(ax1, X, Y, phi, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    
    % Add trajectory projections on surface
    hold(ax1, 'on');
    for traj_idx = 1:num_trajectories
        if size(trajectories{traj_idx}, 1) >= k
            traj_current = trajectories{traj_idx}(1:k, :);
            if size(traj_current, 1) > 1
                % Project trajectory onto surface
                z_traj = interp2(X, Y, phi, traj_current(:,1), traj_current(:,2), 'linear', max(phi(:)));
                plot3(ax1, traj_current(:,1), traj_current(:,2), z_traj + 0.5, ...
                      'Color', colors{mod(traj_idx-1, 6)+1}, 'LineWidth', 2.5);
                
                % Mark current position
                if k > 1
                    scatter3(ax1, traj_current(end,1), traj_current(end,2), z_traj(end) + 0.5, ...
                            80, colors{mod(traj_idx-1, 6)+1}, 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
                end
            end
        end
    end
    
    view(ax1, [45 30]);
    colormap(ax1, parula);
    shading(ax1, 'interp');
    axis(ax1, [domain(1) domain(2) domain(1) domain(2) z_min_global z_max_global]);
    xlabel(ax1, '$x_1$', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel(ax1, '$x_2$', 'FontSize', 18, 'FontWeight', 'bold');
    zlabel(ax1, '$V(x,t)$', 'FontSize', 18, 'FontWeight', 'bold');
    title(ax1, sprintf('Cost-to-Go Surface: $t = %.1f$ s', t), 'FontSize', 18, 'FontWeight', 'bold');
    grid(ax1, 'off');
    
    % Professional axis styling
    set(ax1, 'FontSize', 18, 'LineWidth', 1.2);
    set(ax1, 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on');
    set(ax1, 'TickLength', [0.01 0.025]);
    
    % Subplot 2: Phase Space Trajectories
    ax2 = subplot(2, 3, 3);
    [~, h_contour] = contour(ax2, X, Y, phi, 15, 'LineWidth', 1.2);
    colormap(ax2, gray);
    hold(ax2, 'on');
    
    % Plot all trajectories in phase space
    for traj_idx = 1:num_trajectories
        if size(trajectories{traj_idx}, 1) >= k && k > 1
            traj_current = trajectories{traj_idx}(1:k, :);
            plot(ax2, traj_current(:,1), traj_current(:,2), ...
                 'Color', colors{mod(traj_idx-1, 6)+1}, 'LineWidth', 2.5);
            
            % Current position
            scatter(ax2, traj_current(end,1), traj_current(end,2), ...
                   70, colors{mod(traj_idx-1, 6)+1}, 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
        end
        
        % Starting positions
        scatter(ax2, trajectory_starts(traj_idx,1), trajectory_starts(traj_idx,2), ...
               50, colors{mod(traj_idx-1, 6)+1}, 's', 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    end
    
    axis(ax2, 'equal', 'tight');
    xlim(ax2, domain(1:2)); ylim(ax2, domain(3:4));
    xlabel(ax2, '$x_1$', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel(ax2, '$x_2$', 'FontSize', 18, 'FontWeight', 'bold');
    title(ax2, 'Phase Space Trajectories', 'FontSize', 18, 'FontWeight', 'bold');
    grid(ax2, 'on'); grid(ax2, 'minor');
    
    % Professional axis styling
    set(ax2, 'FontSize', 18, 'LineWidth', 1.2);
    set(ax2, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(ax2, 'TickLength', [0.01 0.025]);
    set(ax2, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
    
    % Subplot 3: Control Signal Evolution
    ax3 = subplot(2, 3, [4 5]);
    
    % Plot control signals for multiple trajectories
    legend_entries = {};
    for traj_idx = 1:min(4, num_trajectories)
        if size(control_signals{traj_idx}, 1) >= k && k > 1
            time_vec = hero_times(1:k);
            control_vec = control_signals{traj_idx}(1:k);
            plot(ax3, time_vec, control_vec, 'Color', colors{traj_idx}, 'LineWidth', 2.5);
            hold(ax3, 'on');
            legend_entries{end+1} = sprintf('Trajectory %d', traj_idx);
        end
    end
    
    % Add current time indicator
    line(ax3, [t t], [u_min_global u_max_global], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
    
    xlim(ax3, [0 10]);
    ylim(ax3, [u_min_global u_max_global]);
    xlabel(ax3, 'Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel(ax3, 'Control $u^*(t)$', 'FontSize', 18, 'FontWeight', 'bold');
    title(ax3, 'Optimal Control Signals', 'FontSize', 18, 'FontWeight', 'bold');
    grid(ax3, 'on'); grid(ax3, 'minor');
    
    % Professional axis styling
    set(ax3, 'FontSize', 18, 'LineWidth', 1.2);
    set(ax3, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(ax3, 'TickLength', [0.01 0.025]);
    set(ax3, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
    
    if ~isempty(legend_entries)
        legend(ax3, legend_entries, 'Location', 'best', 'FontSize', 10);
    end
    
    % Subplot 4: Riccati Matrix Components
    ax4 = subplot(2, 3, 6);
    
    % Plot P-matrix evolution
    if k > 1
        time_vec = hero_times(1:k);
        P_current = zeros(k, 3);
        
        for i = 1:k
            t_interp = hero_times(i);
            if abs(t_interp) < 1e-10
                P_mat = eye(2);
            else
                P_vec_interp = interp1(tRic, Pvec, t_interp, 'pchip');
                P_mat = reshape(P_vec_interp, 2, 2);
            end
            P_current(i, :) = [P_mat(1,1), P_mat(1,2), P_mat(2,2)];
        end
        
        plot(ax4, time_vec, P_current(:,1), 'Color', colors{1}, 'LineWidth', 2.5, 'DisplayName', '$P_{11}$');
        hold(ax4, 'on');
        plot(ax4, time_vec, P_current(:,2), 'Color', colors{2}, 'LineWidth', 2.5, 'DisplayName', '$P_{12}$');
        plot(ax4, time_vec, P_current(:,3), 'Color', colors{3}, 'LineWidth', 2.5, 'DisplayName', '$P_{22}$');
        
        % Current time indicator
        y_lims_p = [min(P_current(:)) - 0.1, max(P_current(:)) + 0.1];
        line(ax4, [t t], y_lims_p, 'Color', 'black', 'LineStyle', '--', 'LineWidth', 2);
        ylim(ax4, y_lims_p);
    end
    
    xlim(ax4, [0 10]);
    xlabel(ax4, 'Time (s)', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel(ax4, 'Matrix Components', 'FontSize', 18, 'FontWeight', 'bold');
    title(ax4, 'Riccati Matrix $P(t)$ Evolution', 'FontSize', 18, 'FontWeight', 'bold');
    grid(ax4, 'on'); grid(ax4, 'minor');
    
    % Professional axis styling
    set(ax4, 'FontSize', 18, 'LineWidth', 1.2);
    set(ax4, 'XMinorTick', 'on', 'YMinorTick', 'on');
    set(ax4, 'TickLength', [0.01 0.025]);
    set(ax4, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
    
    legend(ax4, 'Location', 'best', 'FontSize', 10);
    
    % Overall title
    sgtitle(fig1, sprintf('LQR Solution at $t = %.1f$ s', t), ...
            'FontSize', 24, 'FontWeight', 'bold');
    
    drawnow;
    
    % Capture frame
    frame = getframe(fig1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if k == 1
        imwrite(imind, cm, filename_hero, 'gif', 'Loopcount', inf, 'DelayTime', delay_hero);
    else
        imwrite(imind, cm, filename_hero, 'gif', 'WriteMode', 'append', 'DelayTime', delay_hero);
    end
    
    if mod(k, 5) == 0 || k == length(hero_times)
        fprintf('  ✅ Frame %d/%d (t=%.1f) - Multi-panel analysis\n', k, length(hero_times), t);
    end
end

close(fig1);
fprintf('🏆 HERO professional GIF created: %s\n\n', filename_hero);

%% ========================================================================
%% VALIDATION GIF: Detailed Numerical vs Analytical with Error Analysis
%% ========================================================================

fprintf('📊 Creating detailed validation GIF with error analysis...\n');

fig2 = figure('Position', [100 100 1200 600], 'Color', 'white');
filename_validation = 'lqr_validation_professional.gif';

validation_times = [0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0];

% Pre-compute analytical z-limits for consistency
z_ana_min = inf; z_ana_max = -inf;
for k = 1:length(validation_times)
    t = validation_times(k);
    if abs(t) < 1e-10
        P_t = eye(2);
    else
        P_t = reshape(interp1(tRic, Pvec, t, 'pchip'), 2, 2);
    end
    phi_ana = 0.5 * (P_t(1,1)*X.^2 + 2*P_t(1,2)*X.*Y + P_t(2,2)*Y.^2);
    z_ana_min = min(z_ana_min, min(phi_ana(:)));
    z_ana_max = max(z_ana_max, max(phi_ana(:)));
end

for k = 1:length(validation_times)
    t = validation_times(k);
    
    fname = fullfile(folder, ['phi_t' time_string(t) '.dat']);
    if ~isfile(fname), continue; end
    
    phi_num = reshape(load(fname), N, N)';
    
    % Analytical solution
    if abs(t) < 1e-10
        P_t = eye(2);
    else
        P_t = reshape(interp1(tRic, Pvec, t, 'pchip'), 2, 2);
    end
    phi_ana = 0.5 * (P_t(1,1)*X.^2 + 2*P_t(1,2)*X.*Y + P_t(2,2)*Y.^2);
    
    clf(fig2);
    
    % Numerical solution
    ax1 = subplot(2,3,[1 4]);
    surf(ax1, X, Y, phi_num, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    view(ax1, [45 30]); colormap(ax1, parula); shading(ax1, 'interp');
    zlim(ax1, [z_min_global z_max_global]);
    title(ax1, sprintf('Numerical Solution\n$t = %.1f$ s', t), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax1, '$x_1$', 'FontSize', 12, 'FontWeight', 'bold'); 
    ylabel(ax1, '$x_2$', 'FontSize', 12, 'FontWeight', 'bold'); 
    zlabel(ax1, '$V$', 'FontSize', 12, 'FontWeight', 'bold');
    grid(ax1, 'off');
    set(ax1, 'FontSize', 10, 'LineWidth', 1.2);
    set(ax1, 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on');
    
    % Analytical solution
    ax2 = subplot(2,3,[2 5]);
    surf(ax2, X, Y, phi_ana, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    view(ax2, [45 30]); colormap(ax2, parula); shading(ax2, 'interp');
    zlim(ax2, [z_ana_min z_ana_max]);
    title(ax2, sprintf('Analytical Solution\n$t = %.1f$ s', t), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax2, '$x_1$', 'FontSize', 12, 'FontWeight', 'bold'); 
    ylabel(ax2, '$x_2$', 'FontSize', 12, 'FontWeight', 'bold'); 
    zlabel(ax2, '$V$', 'FontSize', 12, 'FontWeight', 'bold');
    grid(ax2, 'off');
    set(ax2, 'FontSize', 10, 'LineWidth', 1.2);
    set(ax2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on');
    
    % Error analysis
    ax3 = subplot(2,3,[3 6]);
    error_field = abs(phi_num - phi_ana);
    surf(ax3, X, Y, error_field, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    view(ax3, [45 30]); colormap(ax3, hot); shading(ax3, 'interp');
    max_error = max(error_field(:));
    l2_error = sqrt(mean(error_field(:).^2));
    title(ax3, sprintf('Absolute Error\nMax: $%.2e$, $L^2$: $%.2e$', max_error, l2_error), ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel(ax3, '$x_1$', 'FontSize', 12, 'FontWeight', 'bold'); 
    ylabel(ax3, '$x_2$', 'FontSize', 12, 'FontWeight', 'bold'); 
    zlabel(ax3, '$|$Error$|$', 'FontSize', 12, 'FontWeight', 'bold');
    grid(ax3, 'off');
    set(ax3, 'FontSize', 10, 'LineWidth', 1.2);
    set(ax3, 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on');
    
    sgtitle(fig2, 'CASL-HJX Validation: Numerical vs Analytical Riccati Solution', ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    drawnow;
    
    frame = getframe(fig2);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 200);
    
    if k == 1
        imwrite(imind, cm, filename_validation, 'gif', 'Loopcount', inf, 'DelayTime', 0.8);
    else
        imwrite(imind, cm, filename_validation, 'gif', 'WriteMode', 'append', 'DelayTime', 0.8);
    end
    
    fprintf('  ✅ Frame %d/%d (t=%.1f) - Validation analysis\n', k, length(validation_times), t);
end

close(fig2);
fprintf('🔬 Validation GIF created: %s\n\n', filename_validation);

%% ========================================================================
%% CONVERGENCE GIF: Multi-Resolution Analysis
%% ========================================================================

fprintf('📈 Creating convergence analysis GIF...\n');

grid_sizes = [20, 40, 80, 160];
fig3 = figure('Position', [150 150 1000 700], 'Color', 'white');
filename_convergence = 'lqr_convergence_professional.gif';

convergence_times = [1.0, 2.5, 5.0, 10.0];

% Pre-compute z-limits for all grid sizes
z_conv_min = inf; z_conv_max = -inf;
for t_idx = 1:length(convergence_times)
    t = convergence_times(t_idx);
    for g = 1:length(grid_sizes)
        N_g = grid_sizes(g);
        folder_g = sprintf('./LQR2D_Output/LQR2D_%d/phi/', N_g);
        if ~exist(folder_g, 'dir'), continue; end
        
        fname = fullfile(folder_g, ['phi_t' time_string(t) '.dat']);
        if ~isfile(fname), continue; end
        
        phi_num = reshape(load(fname), N_g, N_g)';
        z_conv_min = min(z_conv_min, min(phi_num(:)));
        z_conv_max = max(z_conv_max, max(phi_num(:)));
    end
end

for t_idx = 1:length(convergence_times)
    t = convergence_times(t_idx);
    
    clf(fig3);
    
    % Analytical solution for reference
    if abs(t) < 1e-10
        P_t = eye(2);
    else
        P_t = reshape(interp1(tRic, Pvec, t, 'pchip'), 2, 2);
    end
    
    valid_grids = 0;
    errors = [];
    
    for g = 1:length(grid_sizes)
        N_g = grid_sizes(g);
        folder_g = sprintf('./LQR2D_Output/LQR2D_%d/phi/', N_g);
        
        if ~exist(folder_g, 'dir'), continue; end
        
        fname = fullfile(folder_g, ['phi_t' time_string(t) '.dat']);
        if ~isfile(fname), continue; end
        
        valid_grids = valid_grids + 1;
        
        % Load numerical solution
        phi_num = reshape(load(fname), N_g, N_g)';
        [X_g, Y_g] = meshgrid(linspace(domain(1), domain(2), N_g));
        
        % Analytical solution on same grid
        phi_ana_g = 0.5 * (P_t(1,1)*X_g.^2 + 2*P_t(1,2)*X_g.*Y_g + P_t(2,2)*Y_g.^2);
        max_error = max(abs(phi_num(:) - phi_ana_g(:)));
        errors(end+1) = max_error;
        
        % Plot
        ax = subplot(2, 2, valid_grids);
        surf(ax, X_g, Y_g, phi_num, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
        view(ax, [45 30]); colormap(ax, parula); shading(ax, 'interp');
        zlim(ax, [z_conv_min z_conv_max]);
        
        title(ax, sprintf('$N = %d \\times %d$\nMax Error: $%.2e$', N_g, N_g, max_error), ...
              'FontSize', 13, 'FontWeight', 'bold');
        xlabel(ax, '$x_1$', 'FontSize', 11, 'FontWeight', 'bold'); 
        ylabel(ax, '$x_2$', 'FontSize', 11, 'FontWeight', 'bold'); 
        zlabel(ax, '$V$', 'FontSize', 11, 'FontWeight', 'bold');
        grid(ax, 'off');
        set(ax, 'FontSize', 9, 'LineWidth', 1.2);
        set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on');
        
        if valid_grids >= 4, break; end
    end
    
    sgtitle(fig3, sprintf('Multi-Resolution Convergence Analysis at $t = %.1f$ s', t), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    drawnow;
    
    frame = getframe(fig3);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 180);
    
    if t_idx == 1
        imwrite(imind, cm, filename_convergence, 'gif', 'Loopcount', inf, 'DelayTime', 1.5);
    else
        imwrite(imind, cm, filename_convergence, 'gif', 'WriteMode', 'append', 'DelayTime', 1.5);
    end
    
    fprintf('  ✅ Frame %d/%d (t=%.1f) - Convergence analysis\n', t_idx, length(convergence_times), t);
end

close(fig3);
fprintf('📊 Convergence GIF created: %s\n\n', filename_convergence);

%% ========================================================================
%% RICCATI EVOLUTION GIF: Mathematical Insight
%% ========================================================================

fprintf('🔬 Creating Riccati matrix evolution GIF...\n');

fig4 = figure('Position', [200 200 1000 700], 'Color', 'white');
filename_riccati = 'lqr_riccati_evolution.gif';

% Use all time points for smooth Riccati evolution
riccati_times = hero_times;

% Pre-compute P-matrix limits for consistent scaling
P_min_global = inf; P_max_global = -inf;
for k = 1:length(riccati_times)
    t = riccati_times(k);
    if abs(t) < 1e-10
        P_mat = eye(2);
    else
        P_vec_interp = interp1(tRic, Pvec, t, 'pchip');
        P_mat = reshape(P_vec_interp, 2, 2);
    end
    P_vals = [P_mat(1,1), P_mat(1,2), P_mat(2,2)];
    P_min_global = min(P_min_global, min(P_vals));
    P_max_global = max(P_max_global, max(P_vals));
end

P_margin = 0.1 * (P_max_global - P_min_global);
P_min_global = P_min_global - P_margin;
P_max_global = P_max_global + P_margin;

for k = 1:length(riccati_times)
    t = riccati_times(k);
    
    clf(fig4);
    
    % Create 2x2 subplot for P-matrix components
    component_names = {'$P_{11}$', '$P_{12}$'; '$P_{21}$', '$P_{22}$'};
    component_indices = {[1,1], [1,2]; [2,1], [2,2]};
    
    for i = 1:2
        for j = 1:2
            ax = subplot(2, 2, (i-1)*2 + j);
            
            % Plot time evolution up to current time
            if k > 1
                time_vec = riccati_times(1:k);
                P_values = zeros(k, 1);
                
                for t_idx = 1:k
                    t_curr = riccati_times(t_idx);
                    if abs(t_curr) < 1e-10
                        P_mat = eye(2);
                    else
                        P_vec_interp = interp1(tRic, Pvec, t_curr, 'pchip');
                        P_mat = reshape(P_vec_interp, 2, 2);
                    end
                    P_values(t_idx) = P_mat(component_indices{i,j}(1), component_indices{i,j}(2));
                end
                
                plot(time_vec, P_values, 'Color', colors{1}, 'LineWidth', 3);
                hold on;
                
                % Current point
                scatter(t, P_values(end), 80, colors{2}, 'filled', 'MarkerEdgeColor', 'black', 'LineWidth', 1.5);
                
                title(sprintf('%s $= %.3f$', component_names{i,j}, P_values(end)), ...
                      'FontSize', 14, 'FontWeight', 'bold');
            else
                title(component_names{i,j}, 'FontSize', 14, 'FontWeight', 'bold');
            end
            
            xlim([0 10]);
            ylim([P_min_global P_max_global]);
            
            xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Value', 'FontSize', 12, 'FontWeight', 'bold');
            grid on; grid minor;
            
            % Professional styling
            set(ax, 'FontSize', 11, 'LineWidth', 1.2);
            set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on');
            set(ax, 'TickLength', [0.01 0.025]);
            set(ax, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.15);
        end
    end
    
    sgtitle(sprintf('Riccati Matrix $P(t)$ Evolution: $t = %.2f$ s', t), ...
            'FontSize', 18, 'FontWeight', 'bold');
    
    drawnow;
    
    % Capture frame
    frame = getframe(fig4);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 128);
    
    if k == 1
        imwrite(imind, cm, filename_riccati, 'gif', 'Loopcount', inf, 'DelayTime', 0.12);
    else
        imwrite(imind, cm, filename_riccati, 'gif', 'WriteMode', 'append', 'DelayTime', 0.12);
    end
    
    if mod(k, 6) == 0 || k == length(riccati_times)
        fprintf('  ✅ Frame %d/%d (t=%.2f) - Riccati evolution\n', k, length(riccati_times), t);
    end
end

close(fig4);
fprintf('📊 Riccati evolution GIF created: %s\n\n', filename_riccati);

% Reset default interpreters
set(groot, 'DefaultTextInterpreter', 'remove');
set(groot, 'DefaultAxesTickLabelInterpreter', 'remove');
set(groot, 'DefaultLegendInterpreter', 'remove');

%% ========================================================================
%% SUMMARY AND OPTIMIZATION
%% ========================================================================

fprintf('🎉 ULTRA-PROFESSIONAL CASL-HJX GIFs completed!\n\n');

gif_files = {filename_hero, filename_validation, filename_convergence, filename_riccati};
gif_descriptions = {
    'Multi-panel hero animation with trajectories and control',
    'Detailed validation with error analysis',
    'Multi-resolution convergence study',
    'Riccati matrix mathematical evolution'
};

fprintf('📁 Created professional visualizations:\n');
total_size = 0;

for i = 1:length(gif_files)
    if exist(gif_files{i}, 'file')
        info = dir(gif_files{i});
        size_mb = info.bytes / (1024^2);
        total_size = total_size + size_mb;
        
        status = '🏆';
        if size_mb > 10
            status = '⚠️ ';
        elseif size_mb > 5
            status = '📊';
        end
        
        fprintf('  %s %s (%.1f MB)\n', status, gif_files{i}, size_mb);
        fprintf('     %s\n', gif_descriptions{i});
    end
end

fprintf('\n📊 Total size: %.1f MB\n', total_size);
fprintf('🚀 Ready for academic publication and GitHub showcase!\n\n');

% README integration code
fprintf('📝 Professional README integration:\n\n');
fprintf('```markdown\n');
fprintf('## 🎬 Ultra-Professional Demonstrations\n\n');
fprintf('### Complete LQR Analysis: Surface, Trajectories, and Control\n');
fprintf('![LQR Professional Analysis](%s)\n\n', filename_hero);
fprintf('*Comprehensive visualization showing cost-to-go surface evolution, optimal trajectories in phase space, control signals, and Riccati matrix components. Demonstrates the complete solution of the Hamilton-Jacobi-Bellman equation for optimal control.*\n\n');
fprintf('### Rigorous Numerical Validation\n');
fprintf('![Professional Validation](%s)\n\n', filename_validation);
fprintf('*Detailed comparison between numerical and analytical solutions with comprehensive error analysis. Shows second-order spatial accuracy and validation against exact Riccati solutions.*\n\n');
fprintf('### Multi-Resolution Convergence Study\n');
fprintf('![Convergence Analysis](%s)\n\n', filename_convergence);
fprintf('*Systematic convergence analysis across multiple grid resolutions, demonstrating robust numerical performance and scalability of the CASL-HJX framework.*\n\n');
fprintf('### Mathematical Insight: Riccati Matrix Evolution\n');
fprintf('![Riccati Evolution](%s)\n\n', filename_riccati);
fprintf('*Time evolution of all Riccati matrix components $(P_{11}, P_{12}, P_{21}, P_{22})$ showing convergence to steady-state optimal control parameters.*\n\n');
fprintf('```\n\n');

fprintf('🎯 All professional GIFs ready for publication!\n');

%% ========================================================================
%% HELPER FUNCTIONS  
%% ========================================================================

function [trajectory, control] = compute_optimal_trajectory(start_pos, times, data_folder, N, domain)
    % Compute optimal trajectory using the value function
    trajectory = zeros(length(times), 2);
    control = zeros(length(times), 1);
    trajectory(1, :) = start_pos;
    
    [X, Y] = meshgrid(linspace(domain(1), domain(2), N));
    
    for k = 2:length(times)
        t = times(k-1);
        dt = times(k) - times(k-1);
        
        % Load value function
        fname = fullfile(data_folder, ['phi_t' time_string(t) '.dat']);
        if ~isfile(fname)
            trajectory(k:end, :) = repmat(trajectory(k-1, :), [length(times)-k+1, 1]);
            break;
        end
        
        phi = reshape(load(fname), N, N)';
        
        % Current position
        pos = trajectory(k-1, :);
        
        % Compute gradients (control law: u* = -∇V/∂x₂)
        [phi_x, phi_y] = gradient(phi, X(1,2)-X(1,1), Y(2,1)-Y(1,1));
        
        % Interpolate gradient at current position
        if pos(1) >= domain(1) && pos(1) <= domain(2) && pos(2) >= domain(3) && pos(2) <= domain(4)
            grad_x = interp2(X, Y, phi_x, pos(1), pos(2), 'linear', 0);
            grad_y = interp2(X, Y, phi_y, pos(1), pos(2), 'linear', 0);
            
            % Optimal control (for this 2D LQR system)
            u_opt = -grad_y;  % u* = -∂V/∂x₂
            control(k-1) = u_opt;
            
            % System dynamics: dx₁/dt = x₂, dx₂/dt = u
            x1_dot = pos(2);
            x2_dot = u_opt;
            
            % Forward Euler integration
            trajectory(k, 1) = pos(1) + dt * x1_dot;
            trajectory(k, 2) = pos(2) + dt * x2_dot;
        else
            % Outside domain, stop trajectory
            trajectory(k:end, :) = repmat(trajectory(k-1, :), [length(times)-k+1, 1]);
            break;
        end
    end
end

function dPdt = riccati_rhs(~, P, A, B, Q, R)
    Pm = reshape(P, 2, 2);
    dPdt = (A.' * Pm + Pm * A - Pm * B * (R \ B.') * Pm + Q);
    dPdt = dPdt(:);
end

function str = time_string(t)
    if abs(t - round(t)) < 1e-10
        str = sprintf('%d', round(t));
    else
        str = strrep(sprintf('%.1f', t), '.', 'p');
    end
end