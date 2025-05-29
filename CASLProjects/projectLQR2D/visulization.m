%% GitHub-Optimized LQR2D GIF Creator
%  Creates smaller, web-optimized GIFs perfect for GitHub README

clear; close all; clc;

fprintf('🚀 Creating GitHub-optimized GIFs for LQR2D solver...\n\n');

% Configuration for GitHub optimization
domain = [-4 4 -4 4];
N = 80;  % Reduced resolution for smaller file size
folder = sprintf('./LQR2D_Output/LQR2D_%d/phi/', N);

% Fallback to 160 if 80 doesn't exist
if ~exist(folder, 'dir')
    N = 160;
    folder = sprintf('./LQR2D_Output/LQR2D_%d/phi/', N);
end

if ~exist(folder, 'dir')
    error('No LQR2D data found. Please run the solver first with N=80 or N=160.');
end

fprintf('📊 Using data from: %s\n', folder);

% Optimized time points (fewer frames = smaller file)
gif_times = [0 0.2 0.5 1.0 2.0 3.0 5.0 8.0 10.0];
[X, Y] = meshgrid(linspace(domain(1), domain(2), N));

% Analytical reference
A = [0 1; 0 0]; B = [0; 1]; Q = eye(2); R = 1;
eye_2 = eye(2); 
[tRic, Pvec] = ode45(@(t,P) riccati_rhs(t,P,A,B,Q,R), [0 10], eye_2(:));

%% ========================================================================
%% HERO GIF: Surface Evolution (Main README showcase)
%% ========================================================================

fprintf('🎬 Creating HERO surface evolution GIF...\n');

% Optimized figure size for web
fig = figure('Position', [100 100 600 450], 'Color', 'white', 'MenuBar', 'none');
set(fig, 'InvertHardcopy', 'off');

ax = axes('Position', [0.1 0.15 0.8 0.75]);

filename = 'lqr_solver_demo.gif';
delay = 0.5;  % Longer delay for readability

for k = 1:length(gif_times)
    t = gif_times(k);
    
    fname = fullfile(folder, ['phi_t' time_string(t) '.dat']);
    if ~isfile(fname), continue; end
    
    phi = reshape(load(fname), N, N)';
    
    % Clean, professional 3D plot
    cla(ax);
    surf(ax, X, Y, phi, 'EdgeColor', 'none', 'FaceAlpha', 0.95);
    
    view(ax, [45 30]);
    colormap(ax, parula);
    shading(ax, 'interp');
    
    % Consistent limits for smooth animation
    axis(ax, [domain(1) domain(2) domain(1) domain(2) 0 20]);
    
    % Clean labeling
    xlabel(ax, 'x₁', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, 'x₂', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel(ax, 'Cost V(x,t)', 'FontSize', 12, 'FontWeight', 'bold');
    
    title(ax, sprintf('LQR Cost-to-Go: t = %.1f', t), ...
          'FontSize', 14, 'FontWeight', 'bold');
    
    % Subtle grid
    grid(ax, 'on');
    set(ax, 'GridAlpha', 0.3, 'FontSize', 10);
    
    drawnow;
    
    % Capture and compress
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 128);  % Reduced colors for smaller size
    
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delay);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
    end
    
    fprintf('  ✓ Frame %d/%d (t=%.1f)\n', k, length(gif_times), t);
end

close(fig);
fprintf('✅ Main demo GIF: %s\n\n', filename);

%% ========================================================================
%% VALIDATION GIF: Numerical vs Analytical (Compact comparison)
%% ========================================================================

fprintf('📈 Creating validation comparison GIF...\n');

fig = figure('Position', [150 150 800 400], 'Color', 'white');

filename_val = 'lqr_validation_compact.gif';
comparison_times = [0 1.0 3.0 5.0 10.0];  % Key moments only

for k = 1:length(comparison_times)
    t = comparison_times(k);
    
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
    
    % Side-by-side comparison
    clf;
    
    subplot(1,2,1);
    contourf(X, Y, phi_num, 20, 'LineColor', 'none');
    colormap(parula); axis equal tight;
    title(sprintf('Numerical t=%.1f', t), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('x₁'); ylabel('x₂');
    
    subplot(1,2,2);
    contourf(X, Y, phi_ana, 20, 'LineColor', 'none');
    colormap(parula); axis equal tight;
    title(sprintf('Analytical t=%.1f', t), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('x₁'); ylabel('x₂');
    
    % Error metric
    max_error = max(abs(phi_num(:) - phi_ana(:)));
    sgtitle(sprintf('CASL-HJX Validation | Max Error: %.2e', max_error), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    drawnow;
    
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 128);
    
    if k == 1
        imwrite(imind, cm, filename_val, 'gif', 'Loopcount', inf, 'DelayTime', 1.0);
    else
        imwrite(imind, cm, filename_val, 'gif', 'WriteMode', 'append', 'DelayTime', 1.0);
    end
    
    fprintf('  ✓ Frame %d/%d (t=%.1f)\n', k, length(comparison_times), t);
end

close(fig);
fprintf('✅ Validation GIF: %s\n\n', filename_val);

%% ========================================================================
%% CAPABILITIES GIF: Multi-solver showcase (if other projects exist)
%% ========================================================================

fprintf('🔧 Creating capabilities showcase...\n');

% Check for other project outputs
projects = {'Advection', 'Burgers', 'AdvectionDiffusion', 'Diffusion'};
capabilities_data = {};

for i = 1:length(projects)
    proj_folder = sprintf('./CASLProjects/project%s/__Output/', projects{i});
    if exist(proj_folder, 'dir')
        % Find a result file
        files = dir(fullfile(proj_folder, '*.dat'));
        if ~isempty(files)
            capabilities_data{end+1} = {projects{i}, proj_folder, files(end).name};
        end
    end
end

if length(capabilities_data) >= 2
    fig = figure('Position', [200 200 800 600], 'Color', 'white');
    
    filename_cap = 'casl_capabilities.gif';
    
    % Create capability showcase frames
    for cap_idx = 1:length(capabilities_data)
        clf;
        
        proj_name = capabilities_data{cap_idx}{1};
        proj_folder = capabilities_data{cap_idx}{2};
        proj_file = capabilities_data{cap_idx}{3};
        
        % Load and display the data
        try
            full_path = fullfile(proj_folder, proj_file);
            data = load(full_path);
            
            % Determine grid size
            n_points = length(data);
            n_side = round(sqrt(n_points));
            
            if n_side^2 == n_points
                % 2D data
                phi = reshape(data, n_side, n_side)';
                [X_proj, Y_proj] = meshgrid(linspace(-4, 4, n_side));
                
                surf(X_proj, Y_proj, phi, 'EdgeColor', 'none');
                view([45 30]); colormap(parula); shading interp;
                
                title(sprintf('CASL-HJX: %s Solver', proj_name), ...
                      'FontSize', 16, 'FontWeight', 'bold');
                xlabel('x'); ylabel('y'); zlabel('φ');
                axis tight; grid off;
                
                % Add description
                descriptions = containers.Map();
                descriptions('Advection') = 'Linear transport equation';
                descriptions('Burgers') = 'Nonlinear Burgers equation';
                descriptions('AdvectionDiffusion') = 'Coupled transport-diffusion';
                descriptions('Diffusion') = 'Heat/diffusion equation';
                descriptions('LQR2D') = 'Optimal control (HJB equation)';
                
                if descriptions.isKey(proj_name)
                    text(0.02, 0.98, descriptions(proj_name), ...
                         'Units', 'normalized', 'FontSize', 12, ...
                         'BackgroundColor', 'white', 'EdgeColor', 'black');
                end
            end
        catch
            % Skip if data can't be loaded
            continue;
        end
        
        drawnow;
        
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 128);
        
        if cap_idx == 1
            imwrite(imind, cm, filename_cap, 'gif', 'Loopcount', inf, 'DelayTime', 2.0);
        else
            imwrite(imind, cm, filename_cap, 'gif', 'WriteMode', 'append', 'DelayTime', 2.0);
        end
        
        fprintf('  ✓ Capability %d/%d (%s)\n', cap_idx, length(capabilities_data), proj_name);
    end
    
    close(fig);
    fprintf('✅ Capabilities GIF: %s\n\n', filename_cap);
else
    fprintf('⚠️  Skipping capabilities GIF (need multiple solvers)\n\n');
end

%% ========================================================================
%% FILE SIZE CHECK AND SUMMARY
%% ========================================================================

fprintf('📁 Final file summary:\n');

gif_files = {'lqr_solver_demo.gif', 'lqr_validation_compact.gif'};
if exist('casl_capabilities.gif', 'file')
    gif_files{end+1} = 'casl_capabilities.gif';
end

total_size = 0;
for i = 1:length(gif_files)
    if exist(gif_files{i}, 'file')
        info = dir(gif_files{i});
        size_mb = info.bytes / (1024^2);
        total_size = total_size + size_mb;
        
        status = '✅';
        note = '';
        if size_mb > 3
            status = '⚠️ ';
            note = ' (large for GitHub)';
        end
        
        fprintf('  %s %s: %.1f MB%s\n', status, gif_files{i}, size_mb, note);
    end
end

fprintf('  📊 Total size: %.1f MB\n\n', total_size);

% README template
fprintf('🚀 Add to your README.md:\n\n');
fprintf('```markdown\n');
fprintf('## 🎬 Live Demonstrations\n\n');
fprintf('### Hamilton-Jacobi-Bellman Solver in Action\n');
fprintf('![LQR Solver Demo](lqr_solver_demo.gif)\n\n');
fprintf('*Real-time solution of optimal control problem showing cost-to-go function V(x,t) evolution*\n\n');
fprintf('### Numerical Validation\n');
fprintf('![Validation](lqr_validation_compact.gif)\n\n');
fprintf('*Comparison between numerical solution and analytical Riccati solution*\n\n');
if exist('casl_capabilities.gif', 'file')
    fprintf('### Multi-Solver Capabilities\n');
    fprintf('![CASL Capabilities](casl_capabilities.gif)\n\n');
    fprintf('*CASL-HJX framework solving multiple PDE types*\n\n');
end
fprintf('```\n');

fprintf('🎉 GitHub-optimized GIFs ready for upload!\n');

%% Helper function
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