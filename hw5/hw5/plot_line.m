
clear; clf;

%% Setup Parameters
timeSteps   = [3831, 5108, 6385];               % Time steps to process

% Updated markers and line styles
markers     = {'d', 'p', 'h'};                  % Diamond, pentagram, hexagram
colors      = {[0.1, 0.5, 0.8], [0.9, 0.4, 0.1], [0.3, 0.7, 0.3]};  % Custom RGB colors
lineStyles  = {'-', '-.', ':'};                % Solid, dash-dot, dotted
markerSpacings = [8, 12, 6];                   % Vary marker spacing for each line

figure('Color', 'w'); hold on;

%% Loop Over Time Steps
for i = 1:length(timeSteps)
    tid       = timeSteps(i);
    color     = colors{i};
    marker    = markers{i};
    lineStyle = lineStyles{i};
    spacing   = markerSpacings(i);

    %% --- Load Serial Data ---
    fileSerial = sprintf('T_x_y_%06d.dat', tid);
    dataSerial = readmatrix(fileSerial);
    sortedSerial = sortrows(dataSerial, [2, 1]);

    xSerial = unique(sortedSerial(:,1));
    ySerial = unique(sortedSerial(:,2));
    nx = length(xSerial);
    ny = length(ySerial);
    TSerial = reshape(sortedSerial(:,3), [ny, nx]);

    %% --- Load Parallel Data ---
    filePattern = sprintf('T_x_y_%06d_*.dat', tid);
    parallelFiles = dir(filePattern);

    combined = [];
    for f = 1:length(parallelFiles)
        filePath = fullfile(parallelFiles(f).folder, parallelFiles(f).name);
        combined = [combined; readmatrix(filePath)];
    end

    % Average duplicate (x,y) points
    [uniqueXY, ~, idx] = unique(combined(:,1:2), 'rows');
    Tavg = accumarray(idx, combined(:,3), [], @mean);
    sortedPar = sortrows([uniqueXY, Tavg], [2, 1]);

    xPar = unique(sortedPar(:,1));
    yPar = unique(sortedPar(:,2));

    % Validate grid consistency
    if numel(sortedPar(:,3)) ~= length(xPar) * length(yPar)
        error('Mismatch in reshaped parallel grid size at time = %d.', tid);
    end

    TPar = reshape(sortedPar(:,3), [length(yPar), length(xPar)]);

    %% --- Midline Interpolation at y = 0.5 ---
    yMid = 0.5;
    for j = 1:ny-1
        if ySerial(j) <= yMid && ySerial(j+1) >= yMid
            j1 = j;
            j2 = j + 1;
            break;
        end
    end
    w = (yMid - ySerial(j1)) / (ySerial(j2) - ySerial(j1));
    TmidSerial = (1 - w) * TSerial(j1, :) + w * TSerial(j2, :);
    TmidPar    = (1 - w) * TPar(j1, :) + w * TPar(j2, :);

    markerIndices = 1:spacing:length(xSerial);

    %% --- Plot Serial Data ---
    plot(xSerial, TmidSerial, lineStyle, ...
        'Color', color, 'LineWidth', 2, ...
        'Marker', marker, 'MarkerIndices', markerIndices, ...
        'MarkerSize', 7, ...
        'DisplayName', sprintf('Serial t=%d', tid));

    %% --- Plot Parallel Data ---
    plot(xSerial, TmidPar, '--', ...
        'Color', color * 0.7, 'LineWidth', 2, ...
        'Marker', marker, 'MarkerIndices', markerIndices + 1, ...
        'MarkerSize', 7, ...
        'DisplayName', sprintf('Parallel t=%d', tid));
end

%% Finalize Plot
xlabel('x', 'FontWeight', 'bold');
ylabel('Temperature (T)', 'FontWeight', 'bold');
title('Mid-y Temperature Profiles for Serial and Parallel Runs', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest');
grid on;
grid minor;
xlim([-0.05 1.05]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');

% Save high-quality figure
print('styled_midy_T_profiles', '-dpng', '-r300');

