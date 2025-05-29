%% Simple CASL-HJX GIF Creator
% Run this from the CASL-HJX root directory

clear; close all; clc;

fprintf('🎬 Creating CASL-HJX Visualizations...\n');

%% Create Hero Animation
create_hero_gif();

%% Create Capabilities Demo
create_capabilities_gif();

fprintf('\n✅ Visualizations complete!\n');
fprintf('Files created: casl_hjx_hero.gif, casl_capabilities.gif\n');

function create_hero_gif()
    fprintf('Creating hero animation...\n');
    
    % Parameters
    times = [0, 1, 2, 4, 6, 10];
    domain = [-3, 3, -3, 3];
    N = 80;
    [X, Y] = meshgrid(linspace(domain(1), domain(2), N));
    
    % Setup figure
    fig = figure('Position', [100 100 800 600], 'Color', 'k');
    filename = 'casl_hjx_hero.gif';
    
    % Analytical LQR parameters
    A = [0 1; 0 0]; B = [0; 1]; Q = eye(2); R = 1;
    
    for k = 1:length(times)
        t = times(k);
        
        % Simple analytical approximation for demo
        if t == 0
            P_t = eye(2);
        else
            lambda = 1.2;
            p11 = 1 + lambda*t + 0.3*t^2;
            p12 = 0.4*t + 0.2*t^2;
            p22 = 1 + 0.4*t;
            P_t = [p11, p12; p12, p22];
        end
        
        % Generate cost-to-go surface
        phi = zeros(N, N);
        for i = 1:N
            for j = 1:N
                x_vec = [X(i,j); Y(i,j)];
                phi(i,j) = 0.5 * x_vec' * P_t * x_vec;
            end
        end
        
        % Create visualization
        clf;
        surf(X, Y, phi, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
        view([45 35]); colormap(parula); shading interp;
        
        % Styling
        title(sprintf('CASL-HJX: Hamilton-Jacobi-Bellman Solver (t = %.1f)', t), ...
              'FontSize', 18, 'FontWeight', 'bold', 'Color', 'w');
        xlabel('State x₁', 'FontSize', 14, 'Color', 'w');
        ylabel('State x₂', 'FontSize', 14, 'Color', 'w');
        zlabel('Cost-to-Go V(x,t)', 'FontSize', 14, 'Color', 'w');
        
        xlim(domain(1:2)); ylim(domain(3:4)); zlim([0 15]);
        set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
        
        % Capture frame
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
    fprintf('   ✓ Hero animation saved: %s\n', filename);
end

function create_capabilities_gif()
    fprintf('Creating capabilities demo...\n');
    
    fig = figure('Position', [100 100 1000 600], 'Color', 'w');
    filename = 'casl_capabilities.gif';
    
    % Different PDE examples
    examples = {
        'Hamilton-Jacobi-Bellman', @demo_hjb;
        'Advection Equation', @demo_advection;
        'Diffusion Equation', @demo_diffusion;
        'Level-Set Method', @demo_levelset
    };
    
    for k = 1:length(examples)
        clf;
        
        [X, Y, Z] = examples{k,2}();
        
        % Main surface plot
        subplot(1,2,1);
        surf(X, Y, Z, 'EdgeColor', 'none');
        view([45 30]); colormap(parula); shading interp;
        title(examples{k,1}, 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('x'); ylabel('y'); zlabel('φ');
        
        % Contour plot
        subplot(1,2,2);
        contourf(X, Y, Z, 20, 'LineColor', 'none');
        title('Contour View', 'FontSize', 12);
        colorbar; axis equal tight;
        
        sgtitle(sprintf('CASL-HJX: %s', examples{k,1}), 'FontSize', 16, 'FontWeight', 'bold');
        
        drawnow;
        frame = getframe(gcf);
        [imind, cm] = rgb2ind(frame.cdata, 256);
        
        if k == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 2.0);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 2.0);
        end
    end
    
    close(fig);
    fprintf('   ✓ Capabilities demo saved: %s\n', filename);
end

% Helper functions for different equation types
function [X, Y, Z] = demo_hjb()
    [X, Y] = meshgrid(linspace(-2, 2, 60));
    Z = 0.5*(X.^2 + Y.^2) + 0.1*sin(3*X).*cos(3*Y);
end

function [X, Y, Z] = demo_advection()
    [X, Y] = meshgrid(linspace(-2, 2, 60));
    Z = exp(-0.5*((X-0.5).^2 + (Y-0.5).^2));
end

function [X, Y, Z] = demo_diffusion()
    [X, Y] = meshgrid(linspace(-2, 2, 60));
    Z = exp(-(X.^2 + Y.^2)/0.8);
end

function [X, Y, Z] = demo_levelset()
    [X, Y] = meshgrid(linspace(-2, 2, 60));
    Z = tanh(2*(sqrt(X.^2 + Y.^2) - 1));
end