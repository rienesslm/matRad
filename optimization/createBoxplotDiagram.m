function createBoxplotDiagram(boxplotData)
% Input: boxplotData - a matrix where each row represents data for a boxplot
% Each row represents data for a different "Generation" (boxplot).

numGenerations = size(boxplotData, 1);
generations = 1:numGenerations;

% Create the figure and set axis labels
figure;

% Plot the boxplot
boxplot(boxplotData);
xlabel('Generations');
ylabel('Objective values');
title('Optimization Progress over Generations');

% Customize the plot appearance
grid on;
end