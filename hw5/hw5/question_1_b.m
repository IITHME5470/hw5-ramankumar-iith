% Grid dimensions
grid_width = 800;  % nx
grid_height = 800; % ny

% Initialize solution matrices
parallel_temperature = NaN(grid_width, grid_height); 
serial_temperature = NaN(grid_width, grid_height);  


%parallel
parallel_data_files = dir('parallel_*.dat');  % Find all parallel output files

for file = parallel_data_files'
    % Load data from current file
    file_data = load(file.name);
    
    % Convert from 0-based to 1-based indexing
    x_coords = file_data(:,1) + 1; 
    y_coords = file_data(:,2) + 1;
    temperatures = file_data(:,3);
    
    % Store temperatures in the appropriate grid locations
    for entry = 1:length(x_coords)
        parallel_temperature(x_coords(entry), y_coords(entry)) = temperatures(entry);
    end
end

%serial
serial_data = load('serial_000010.dat');

% Convert from 0-based to 1-based indexing
serial_x = serial_data(:,1) + 1;  
serial_y = serial_data(:,2) + 1;
serial_T = serial_data(:,3);

% Store temperatures in the appropriate grid locations
for entry = 1:length(serial_x)
    serial_temperature(serial_x(entry), serial_y(entry)) = serial_T(entry);
end

% =============================================
% Compare solutions
% =============================================
temperature_differences = abs(serial_temperature - parallel_temperature);

% Calculate statistics
max_difference = max(temperature_differences(:));
min_difference = min(temperature_differences(:));
mean_difference = mean(temperature_differences(:), 'omitnan');

% =============================================
% Display comparison results
% =============================================
fprintf('\nSerial vs. Parallel Solution Comparison\n');
fprintf('--------------------------------------\n');
fprintf('Grid dimensions: %d x %d\n', grid_width, grid_height);
fprintf('Maximum absolute difference: %.5e\n', max_difference);
fprintf('Minimum absolute difference: %.5e\n', min_difference);
fprintf('Mean absolute difference:    %.5e\n', mean_difference);
fprintf('\n');

% Check against machine precision
if max_difference < eps
    fprintf('All differences are within machine precision (eps = %.2e)\n', eps);
else
    fprintf('Differences exceed machine precision! Maximum difference: %.2e\n', max_difference);

end

% Additional quality checks
if any(isnan(temperature_differences))
    missing_points = sum(isnan(temperature_differences(:)));
    fprintf('\nWarning: %d grid points have no data in one or both solutions\n', missing_points);
end