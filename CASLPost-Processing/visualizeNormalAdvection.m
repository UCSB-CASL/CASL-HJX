function visualizeNormalAdvection()
% VISUALIZENORMALADVECTION - Comprehensive visualization for Normal Advection Equation
%
% This script creates:
% 1. Zero level set contour evolution (boundary tracking)
% 2. 3D surface plots of the advected property φ(x,y,t)
% 3. Animation of the evolution process
%
% MATHEMATICAL MODEL: ∂φ/∂t + vₙ|∇φ| = 0
% Initial condition: φ(x,y,0) = 1 - x⁴ - y⁴ + 3x²
% Domain: [-2, 2] × [-2, 2]
%
% Author: Generated for CASL-HJX Normal Advection Project
% Date: May 29, 2025

    %% Configuration and Setup
    fprintf('CASL-HJX Normal Advection Visualization\n');
    fprintf('======================================\n');
    
    % Project paths - Auto-detect current location
    currentDir = pwd;
    if contains(currentDir, 'projectNormalAdvection')
        outputDir = '__Output';
        saveDir = 'PostProcessing';
    else
        % Assume we're in CASLHJB2D root or elsewhere
        outputDir = '../CASLProjects/projectNormalAdvection/__Output';
        saveDir = '../CASLProjects/projectNormalAdvection/PostProcessing';
    end
    
    fprintf('Data directory: %s\n', outputDir);
    fprintf('Output directory: %s\n', saveDir);
    
    % Create output directory for plots
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
        fprintf('Created directory: %s\n', saveDir);
    end
    
    %% Domain Parameters (matching C++ simulation)
    xMin = -2.0; xMax = 2.0;
    yMin = -2.0; yMax = 2.0;
    nX = 250; nY = 250;
    
    % Create coordinate arrays
    x = linspace(xMin, xMax, nX);
    y = linspace(yMin, yMax, nY);
    [X, Y] = meshgrid(x, y);
    
    fprintf('Domain: [%.1f, %.1f] × [%.1f, %.1f]\n', xMin, xMax, yMin, yMax);
    fprintf('Grid: %d×%d points\n', nX, nY);
    
    %% Load Data Files
    fprintf('\nLoading simulation data...\n');
    
    % Find all phi_t*.dat files
    dataFiles = dir(fullfile(outputDir, 'phi_t*.dat'));
    if isempty(dataFiles)
        error('No data files found in %s', outputDir);
    end
    
    % Sort files by time (natural sorting)
    fileNames = {dataFiles.name};
    times = extractTimeFromFilename(fileNames);
    [times, sortIdx] = sort(times);
    fileNames = fileNames(sortIdx);
    
    nFiles = length(fileNames);
    fprintf('Found %d data files\n', nFiles);
    fprintf('Time range: t = %.3f to %.3f\n', times(1), times(end));
    
    %% Load and Process Data
    phi_data = cell(nFiles, 1);
    
    fprintf('Loading data files...\n');
    for i = 1:nFiles
        filename = fullfile(outputDir, fileNames{i});
        phi_data{i} = loadPhiData(filename, nX, nY);
        
        if mod(i, 10) == 0 || i == 1 || i == nFiles
            fprintf('  Loaded %d/%d files (t = %.3f)\n', i, nFiles, times(i));
        end
    end
    
    %% Visualization 1: Zero Level Set Evolution
    fprintf('\nCreating zero level set visualization...\n');
    createZeroLevelSetPlot(X, Y, phi_data, times, saveDir);
    
    %% Visualization 2: 3D Surface Evolution
    fprintf('Creating 3D surface visualization...\n');
    create3DSurfacePlot(X, Y, phi_data, times, saveDir);
    
    %% Visualization 3: Combined Animation
    fprintf('Creating evolution animation...\n');
    createEvolutionAnimation(X, Y, phi_data, times, saveDir);
    
    %% Analysis: Level Set Statistics
    fprintf('Performing level set analysis...\n');
    analyzeLevelSetEvolution(X, Y, phi_data, times, saveDir);
    
    fprintf('\nVisualization complete!\n');
    fprintf('Results saved in: %s\n', saveDir);
    
end

%% Helper Functions

function times = extractTimeFromFilename(fileNames)
    % Extract time values from phi_t*.dat filenames
    times = zeros(size(fileNames));
    
    for i = 1:length(fileNames)
        name = fileNames{i};
        % Remove phi_t and .dat
        timeStr = strrep(name, 'phi_t', '');
        timeStr = strrep(timeStr, '.dat', '');
        
        if strcmp(timeStr, '0')
            times(i) = 0;
        else
            % Replace 'p' with '.' for decimal point
            timeStr = strrep(timeStr, 'p', '.');
            times(i) = str2double(timeStr);
        end
    end
end

function phi = loadPhiData(filename, nX, nY)
    % Load phi data from .dat file (assuming space-separated ASCII format)
    try
        data = dlmread(filename);
        
        % Reshape data to grid format
        if numel(data) == nX * nY
            phi = reshape(data, nY, nX);
        else
            % Try reading as matrix directly
            phi = data;
            if size(phi, 1) ~= nY || size(phi, 2) ~= nX
                error('Data dimensions do not match expected grid size');
            end
        end
    catch ME
        fprintf('Error loading %s: %s\n', filename, ME.message);
        phi = zeros(nY, nX);
    end
end

function createZeroLevelSetPlot(X, Y, phi_data, times, saveDir)
    % Create zero level set evolution plot
    
    figure('Position', [100, 100, 1200, 800]);
    
    % Select representative time steps for display
    nDisplay = min(8, length(times));
    displayIdx = round(linspace(1, length(times), nDisplay));
    
    colors = lines(nDisplay);
    
    hold on;
    legendEntries = cell(nDisplay, 1);
    
    for i = 1:nDisplay
        idx = displayIdx(i);
        phi = phi_data{idx};
        
        % Plot zero level set contour
        contour(X, Y, phi, [0, 0], 'LineWidth', 2, 'Color', colors(i, :));
        legendEntries{i} = sprintf('t = %.3f', times(idx));
    end
    
    xlabel('x');
    ylabel('y');
    title('Zero Level Set Evolution: \partial\phi/\partial t + v_n|\nabla\phi| = 0');
    legend(legendEntries, 'Location', 'best');
    grid on;
    axis equal;
    xlim([min(X(:)), max(X(:))]);
    ylim([min(Y(:)), max(Y(:))]);
    
    % Save plot
    saveas(gcf, fullfile(saveDir, 'zero_levelset_evolution.png'));
    saveas(gcf, fullfile(saveDir, 'zero_levelset_evolution.fig'));
    close(gcf);
end

function create3DSurfacePlot(X, Y, phi_data, times, saveDir)
    % Create 3D surface plots at key time instances
    
    % Select key time points
    nDisplay = min(6, length(times));
    displayIdx = round(linspace(1, length(times), nDisplay));
    
    figure('Position', [100, 100, 1400, 900]);
    
    for i = 1:nDisplay
        subplot(2, 3, i);
        idx = displayIdx(i);
        phi = phi_data{idx};
        
        surf(X, Y, phi, 'EdgeColor', 'none');
        shading interp;
        colorbar;
        
        xlabel('x');
        ylabel('y');
        zlabel('\phi(x,y,t)');
        title(sprintf('t = %.3f', times(idx)));
        
        % Set consistent z-limits for comparison
        if i == 1
            zLimits = [min(phi(:)), max(phi(:))];
        end
        zlim(zLimits);
        
        view(45, 30);
    end
    
    sgtitle('3D Surface Evolution: \phi(x,y,t)');
    
    % Save plot
    saveas(gcf, fullfile(saveDir, 'surface_evolution_3D.png'));
    saveas(gcf, fullfile(saveDir, 'surface_evolution_3D.fig'));
    close(gcf);
end

function createEvolutionAnimation(X, Y, phi_data, times, saveDir)
    % Create professional animated visualization with GIF output
    
    % Set up professional formatting
    figure('Position', [100, 100, 1400, 600], 'Color', 'white');
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    
    % Prepare for animation
    nFrames = length(times);
    frameSkip = max(1, floor(nFrames / 50)); % Limit to ~50 frames for smoother GIF
    animFrames = 1:frameSkip:nFrames;
    
    % Create video writer for MP4
    videoFile = fullfile(saveDir, 'normal_advection_evolution.mp4');
    v = VideoWriter(videoFile, 'MPEG-4');
    v.FrameRate = 8;
    open(v);
    
    % Prepare for GIF creation
    gifFile = fullfile(saveDir, 'normal_advection_evolution.gif');
    
    for frameIdx = 1:length(animFrames)
        i = animFrames(frameIdx);
        clf;
        
        % Left subplot: Zero level set contour
        subplot(1, 2, 1);
        phi = phi_data{i};
        
        contourf(X, Y, phi, 15, 'LineStyle', 'none');
        hold on;
        contour(X, Y, phi, [0, 0], 'Color', [0.8, 0.2, 0.2], 'LineWidth', 3);
        
        cb1 = colorbar('Location', 'eastoutside');
        cb1.TickLabelInterpreter = 'latex';
        cb1.Label.String = '$\phi(x,y,t)$';
        cb1.Label.Interpreter = 'latex';
        cb1.Label.FontSize = 12;
        
        xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
        title(sprintf('$t = %.3f$', times(i)), 'FontSize', 12, 'Interpreter', 'latex');
        axis equal;
        xlim([min(X(:)), max(X(:))]);
        ylim([min(Y(:)), max(Y(:))]);
        
        % Professional tick formatting
        xticks(-2:1:2);
        yticks(-2:1:2);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
        set(gca, 'LineWidth', 1, 'Box', 'on');
        grid off; % Remove grid
        
        % Right subplot: 3D surface
        subplot(1, 2, 2);
        surf(X, Y, phi, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
        shading interp;
        
        cb2 = colorbar('Location', 'eastoutside');
        cb2.TickLabelInterpreter = 'latex';
        cb2.Label.String = '$\phi(x,y,t)$';
        cb2.Label.Interpreter = 'latex';
        cb2.Label.FontSize = 12;
        
        xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
        ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
        zlabel('$\phi(x,y,t)$', 'FontSize', 14, 'Interpreter', 'latex');
        title(sprintf('$t = %.3f$', times(i)), 'FontSize', 12, 'Interpreter', 'latex');
        view(45, 30);
        
        % Consistent z-limits
        if frameIdx == 1
            phi_all = cell2mat(phi_data);
            zLimits = [min(phi_all(:)) - 0.5, max(phi_all(:)) + 0.5];
        end
        zlim(zLimits);
        
        % Professional 3D formatting
        xticks(-2:1:2);
        yticks(-2:1:2);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
        set(gca, 'LineWidth', 1, 'Box', 'on');
        grid off; % Remove grid
        
        % Overall title
        sgtitle('Normal Advection: $\frac{\partial\phi}{\partial t} + v_n|\nabla\phi| = 0$', ...
               'FontSize', 16, 'Interpreter', 'latex', 'FontWeight', 'bold');
        
        drawnow;
        
        % Capture frame for video
        frame = getframe(gcf);
        writeVideo(v, frame);
        
        % Capture frame for GIF
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        
        if frameIdx == 1
            imwrite(imind, cm, gifFile, 'gif', 'Loopcount', inf, 'DelayTime', 0.15);
        else
            imwrite(imind, cm, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
        end
        
        % Progress indicator
        if mod(frameIdx, 10) == 0 || frameIdx == length(animFrames)
            fprintf('  Animation progress: %d/%d frames\n', frameIdx, length(animFrames));
        end
    end
    
    close(v);
    close(gcf);
    
    % Reset MATLAB defaults
    set(0, 'DefaultTextInterpreter', 'remove');
    set(0, 'DefaultAxesTickLabelInterpreter', 'remove');
    
    fprintf('Animation saved: %s\n', videoFile);
    fprintf('GIF saved: %s\n', gifFile);
end

function analyzeLevelSetEvolution(X, Y, phi_data, times, saveDir)
    % Analyze level set evolution statistics
    
    nTimes = length(times);
    zeroArea = zeros(nTimes, 1);
    minPhi = zeros(nTimes, 1);
    maxPhi = zeros(nTimes, 1);
    
    dx = X(1, 2) - X(1, 1);
    dy = Y(2, 1) - Y(1, 1);
    cellArea = dx * dy;
    
    for i = 1:nTimes
        phi = phi_data{i};
        
        % Area inside zero level set (φ > 0)
        zeroArea(i) = sum(phi(:) > 0) * cellArea;
        
        % Min/max values
        minPhi(i) = min(phi(:));
        maxPhi(i) = max(phi(:));
    end
    
    % Create analysis plots
    figure('Position', [100, 100, 1200, 400]);
    
    subplot(1, 3, 1);
    plot(times, zeroArea, 'b-', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Area (\phi > 0)');
    title('Zero Level Set Area');
    grid on;
    
    subplot(1, 3, 2);
    plot(times, minPhi, 'r-', 'LineWidth', 2);
    hold on;
    plot(times, maxPhi, 'b-', 'LineWidth', 2);
    xlabel('Time');
    ylabel('\phi values');
    title('φ Range Evolution');
    legend('min(φ)', 'max(φ)', 'Location', 'best');
    grid on;
    
    subplot(1, 3, 3);
    plot(times, maxPhi - minPhi, 'g-', 'LineWidth', 2);
    xlabel('Time');
    ylabel('Range');
    title('φ Range Width');
    grid on;
    
    sgtitle('Level Set Evolution Analysis');
    
    % Save analysis
    saveas(gcf, fullfile(saveDir, 'levelset_analysis.png'));
    saveas(gcf, fullfile(saveDir, 'levelset_analysis.fig'));
    close(gcf);
    
    % Save numerical data
    analysisData = table(times', zeroArea, minPhi, maxPhi, ...
        'VariableNames', {'Time', 'ZeroLevelSetArea', 'MinPhi', 'MaxPhi'});
    writetable(analysisData, fullfile(saveDir, 'analysis_data.csv'));
    
    fprintf('Analysis complete. Zero level set area change: %.2f%%\n', ...
        100 * (zeroArea(end) - zeroArea(1)) / zeroArea(1));
end