%% LQR Professional GIF Creator
% Run this from the projectLQR2D directory

clear; close all; clc;

fprintf('🎬 Creating Professional LQR Visualizations...\n');

%% Configuration
gridSize = 160;  % Use your 160x160 data
domain = [-4 4 -4 4];  % Your domain
folder = sprintf('./LQR2D_Output/LQR2D_%d/phi/', gridSize);

% Check if data exists
if ~isfolder(folder)
    error('Data folder not found: %s\nMake sure you run this from projectLQR2D directory', folder);
end

% Find available time files
files = dir(fullfile(folder, 'phi_t*.dat'));
if isempty(files)
    error('No phi_t*.dat files found in %s', folder);
end

fprintf('Found %d data files\n', length(files));

%% Extract time points and sort
times = [];
filenames = {};
for i = 1:length(files)
    name = files(i).name;
    % Extract time from filename (e.g., phi_t2p5.dat -> 2.5)
    time_str = extractBetween(name, 'phi_t', '.dat');
    if ~isempty(time_str)
        time_str = strrep(time_str{1}, 'p', '.');
        times(end+1) = str2double(time_str);
        filenames{end+1} = name;
    end
end

% Sort by time
[times, idx] = sort(times);
filenames = filenames(idx);

% Select time points for smooth animation (max 10 frames)
if length(times) > 10
    step = floor(length(times) / 10);
    selected_idx = 1:step:length(times);
    if selected_idx(end) ~= length(times)
        selected_idx(end+1) = length(times);
    end
    times = times(selected_idx);
    filenames = filenames(selected_idx);
end

fprintf('Using %d time points: [', length(times));
fprintf('%.1f ', times);
fprintf(']\n');

%% Setup grid
[X, Y] = meshgrid(linspace(domain(1), domain(2), gridSize));

%% Create Hero GIF - 3D Surface Evolution
fprintf('Creating hero surface animation...\n');
create_surface_gif();

%% Create Contour GIF - 2D Evolution  
fprintf('Creating contour evolution...\n');
create_contour_gif();

%% Create Comparison GIF - Numerical vs Analytical
fprintf('Creating validation comparison...\n');
create_comparison_gif();

fprintf('\n✅ All GIFs created successfully!\n');
fprintf('📁 Files created:\n');
fprintf('   • lqr_hero_3d.gif (Main showcase - 3D surface)\n');
fprintf('   • lqr_contour.gif (2D contour evolution)\n');
fprintf('   • lqr_validation.gif (Numerical vs Analytical)\n\n');
fprintf('💡 Use lqr_hero_3d.gif as your main README animation\n');

%% FUNCTION DEFINITIONS

function create_surface_gif()
    fig = figure('Position', [100 100 800 600], 'Color', 'k');
    filename = 'lqr_hero_3d.gif';
    
    % Find z-limits for consistent scaling
    z_max = 0;
    for k = 1:length(times)
        phi_temp = load_phi_data(k);
        z_max = max(z_max, max(phi_temp(:)));
    end
    
    for k = 1:length(times)
        t = times(k);
        phi = load_phi_data(k);
        
        clf;
        h = surf(X, Y, phi, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
        
        % Professional styling
        view([45 35]);
        colormap(parula);
        shading interp;
        
        % Lighting for premium look
        lighting gouraud;
        light('Position', [1 1 1], 'Style', 'local');
        light('Position', [-1 -1 0.5], 'Style', 'local', 'Color', [0.3 0.3 0.4]);
        
        % Modern typography
        title(sprintf('Hamilton-Jacobi-Bellman Solver: t = %.1f', t), ...
              'FontSize', 20, 'FontWeight', 'bold', 'Color', 'w');
        xlabel('State x₁', 'FontSize', 14, 'Color', 'w', 'FontWeight', 'bold');
        ylabel('State x₂', 'FontSize', 14, 'Color', 'w', 'FontWeight', 'bold');
        zlabel('Cost-to-Go V(x,t)', 'FontSize', 14, 'Color', 'w', 'FontWeight', 'bold');
        
        % Consistent limits
        xlim(domain(1:2)); ylim(domain(3:4)); zlim([0 z_max*1.1]);
        set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', ...
                'FontSize', 12, 'LineWidth', 1.5);
        
        % Branding
        text(0.02, 0.95, 'CASL-HJX Framework', 'Units', 'normalized', ...
             'FontSize', 12, 'Color', [0.7 0.7 0.7], 'FontWeight', 'bold');
        
        % Capture frame
        drawnow;
        frame = getframe(gcf);
        [imind, cm] = rgb2ind(frame.cdata, 256);
        
        if k == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, ...
                    'DelayTime', 0.8, 'Disposal', 'restoreBG');
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', ...
                    'DelayTime', 0.8, 'Disposal', 'restoreBG');
        end
    end
    
    close(fig);
    fprintf('   ✓ 3D surface GIF saved: %s\n', filename);
end

function create_contour_gif()
    fig = figure('Position', [100 100 600 500], 'Color', 'w');
    filename = 'lqr_contour.gif';
    
    % Find color limits
    c_max = 0;
    for k = 1:length(times)
        phi_temp = load_phi_data(k);
        c_max = max(c_max, max(phi_temp(:)));
    end
    
    for k = 1:length(times)
        t = times(k);
        phi = load_phi_data(k);
        
        clf;
        contourf(X, Y, phi, 20, 'LineColor', 'none');
        colormap(parula);
        
        title(sprintf('LQR Cost-to-Go Evolution: t = %.1f', t), ...
              'FontSize', 16, 'FontWeight', 'bold');
        xlabel('State x₁', 'FontSize', 12);
        ylabel('State x₂', 'FontSize', 12);
        
        colorbar;
        caxis([0 c_max]);
        axis equal tight;
        set(gca, 'FontSize', 11);
        
        drawnow;
        frame = getframe(gcf);
        [imind, cm] = rgb2ind(frame.cdata, 256);
        
        if k == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.6);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.6);
        end
    end
    
    close(fig);
    fprintf('   ✓ Contour GIF saved: %s\n', filename);
end

function create_comparison_gif()
    fig = figure('Position', [100 100 1200 500], 'Color', 'w');
    filename = 'lqr_validation.gif';
    
    % Analytical Riccati solution
    A = [0 1; 0 0]; B = [0; 1]; Q = eye(2); R = 1;
    eye_2 = eye(2); 
    [t_ric, P_vec] = ode45(@(t,P) riccati_rhs(t,P,A,B,Q,R), [0 max(times)], eye_2(:));
    
    for k = 1:length(times)
        t = times(k);
        phi_num = load_phi_data(k);
        
        % Compute analytical solution
        P_t = interp1(t_ric, P_vec, t, 'pchip');
        P_mat = reshape(P_t, 2, 2);
        phi_ana = zeros(gridSize, gridSize);
        for i = 1:gridSize
            for j = 1:gridSize
                x_vec = [X(i,j); Y(i,j)];
                phi_ana(i,j) = 0.5 * x_vec' * P_mat * x_vec;
            end
        end
        
        clf;
        
        % Numerical (left)
        subplot(1,3,1);
        contourf(X, Y, phi_num, 20, 'LineColor', 'none');
        title('Numerical Solution', 'FontSize', 12, 'FontWeight', 'bold');
        colorbar; axis equal tight; colormap(parula);
        
        % Analytical (center)
        subplot(1,3,2);
        contourf(X, Y, phi_ana, 20, 'LineColor', 'none');
        title('Analytical Solution', 'FontSize', 12, 'FontWeight', 'bold');
        colorbar; axis equal tight; colormap(parula);
        
        % Error (right)
        subplot(1,3,3);
        error = abs(phi_num - phi_ana);
        contourf(X, Y, error, 20, 'LineColor', 'none');
        title('Absolute Error', 'FontSize', 12, 'FontWeight', 'bold');
        colorbar; axis equal tight; colormap(hot);
        
        sgtitle(sprintf('CASL-HJX Validation: t = %.1f', t), ...
                'FontSize', 16, 'FontWeight', 'bold');
        
        drawnow;
        frame = getframe(gcf);
        [imind, cm] = rgb2ind(frame.cdata, 256);
        
        if k == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1.0);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1.0);
        end
    end
    
    close(fig);
    fprintf('   ✓ Validation GIF saved: %s\n', filename);
end

function phi = load_phi_data(k)
    fname = fullfile(folder, filenames{k});
    phi = reshape(load(fname), gridSize, gridSize)';
end

function dPdt = riccati_rhs(~, P, A, B, Q, R)
    Pm = reshape(P, 2, 2);
    dPdt = (A.' * Pm + Pm * A - Pm * B * (R \ B.') * Pm + Q);
    dPdt = dPdt(:);
end