% Define the font to be used
font = 'Arial';

% Create a figure and set default font properties for axes and text
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

hold on; 
box on;

% Define improved color scheme for clarity
colors = [
    0, 0.4470, 0.7410;   % Blue
    0.8500, 0.3250, 0.0980; % Orange
    0.9290, 0.6940, 0.1250; % Yellow
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.6350, 0.0780, 0.1840; % Maroon
];

% Line styles for variety
lineStyles = {'--', '-.', '-', '--', '-.', '-'};

% Plot CDFs with improved colors and styles
plot(sort(positioning_err_distr_median(1,:),'ascend'), linspace(0,1,nbrOfSetups*num_tp_points), ...
    'Color', colors(1,:), 'LineStyle', lineStyles{1}, 'LineWidth', 2);
plot(sort(positioning_err_distr_median(2,:),'ascend'), linspace(0,1,nbrOfSetups*num_tp_points), ...
    'Color', colors(2,:), 'LineStyle', lineStyles{2}, 'LineWidth', 2);
plot(sort(positioning_err_distr_median(3,:),'ascend'), linspace(0,1,nbrOfSetups*num_tp_points), ...
    'Color', colors(3,:), 'LineStyle', lineStyles{3}, 'LineWidth', 2);
plot(sort(positioning_err_distr_median(4,:),'ascend'), linspace(0,1,nbrOfSetups*num_tp_points), ...
    'Color', colors(4,:), 'LineStyle', lineStyles{4}, 'LineWidth', 2);
plot(sort(positioning_err_distr_median(5,:),'ascend'), linspace(0,1,nbrOfSetups*num_tp_points), ...
    'Color', colors(5,:), 'LineStyle', lineStyles{5}, 'LineWidth', 2);
plot(sort(positioning_err_distr_median(6,:),'ascend'), linspace(0,1,nbrOfSetups*num_tp_points), ...
    'Color', colors(6,:), 'LineStyle', lineStyles{6}, 'LineWidth', 2);

% Add legend with tex interpreter
legend({'K=16', 'K=36', 'K=64', 'K=100', 'K=144', 'K=225'}, ...
    'Interpreter', 'tex', 'Location', 'SouthEast', 'FontSize', 19);

% Set axis properties
set(gca, 'FontSize', 20);

% Add labels with tex interpreter
xlabel('Localization Error - Distributed Median (m)', ...
    'Interpreter', 'tex', 'FontSize', 18, 'FontName', font);
ylabel('CDF', ...
    'Interpreter', 'tex', 'FontSize', 22, 'FontName', font);

xlim([0 25]);

% Add grid
grid on;
set(gca, 'GridLineStyle', ':', 'GridColor', 'k', 'GridAlpha', 0.5);

% Force MATLAB to use 'painters' renderer
set(gcf, 'Renderer', 'painters');

% Save the figure in high resolution
print(gcf, 'Fig9_positioning_err_distr_median_cdf', '-dpng', '-r300'); % Save as PNG
