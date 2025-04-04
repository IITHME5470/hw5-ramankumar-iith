clear
tid = 6385; 

% Find all files for the given timestep with any rank
filePattern = sprintf('T_x_y_%06d_*.dat', tid);
files = dir(filePattern);
filePattern
length(files)

% Read and concatenate data from all ranks
combinedData = [];
for k = 1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    rankData = readmatrix(filePath);
    combinedData = [combinedData; rankData];
end

a = sortrows(combinedData, [2, 1]);

% Extract unique x and y coordinates
x = unique(a(:, 1));  
y = unique(a(:, 2));  

% Determine grid dimensions
nx = length(x);
ny = length(y);

T = reshape(a(:,3), [ny, nx]);

figure, clf
contourf(x,y,T','LineColor', 'none')
xlabel('x'), ylabel('y'), title(strcat('t = ',sprintf(' %06d',tid)));
xlim([-0.05 1.05]), ylim([-0.05 1.05]), clim([-0.05 1.05]), colorbar
colormap('jet')
set(gca, 'FontSize', 14)
screen2jpeg(strcat('cont_T_', sprintf('%04d', tid), '.png'))

% Find the true mid-y (y = 0.5) and interpolate
mid_y = 0.5;
% Find the two y-indices that bracket mid_y
for j = 1:ny-1
    if y(j) <= mid_y && y(j+1) >= mid_y
        j1 = j;
        j2 = j+1;
        break;
    end
end
% Linear interpolation weights
w = (mid_y - y(j1)) / (y(j2) - y(j1));
% Interpolate T at mid_y for each x
Tmid = (1 - w) * T(j1, :) + w * T(j2, :);

figure, clf
plot(x, Tmid, '-', 'LineWidth', 2)
xlabel('x'), ylabel('T'), title(strcat('Profile along mid-y at t=',sprintf('%06d',tid)))
xlim([-0.05 1.05])
set(gca, 'FontSize', 14)
screen2jpeg(strcat('line_midy_T_', sprintf('%04d', tid), '.png'))
