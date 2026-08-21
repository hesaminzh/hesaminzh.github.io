% Display four groups of three bars
figure;
y = [[2,2,3]; [2,5,6]; [2,8,9]; [2,11,12]];
bar(y);



% Data
algorithms = {'Algorithm A', 'Algorithm B', 'Algorithm C'};
test_cases = {'Test 1', 'Test 2', 'Test 3', 'Test 4'};
performance = [85, 90, 95, 80; 88, 85, 92, 86; 90, 87, 93, 89];

% Create grouped bar graph
figure
b = bar(performance', 'grouped');
set(gca, 'XTickLabel', test_cases)
legend(algorithms, 'Location', 'northwest')
xlabel('Test Cases')
ylabel('Performance (%)')
title('Algorithm Performance Comparison')
