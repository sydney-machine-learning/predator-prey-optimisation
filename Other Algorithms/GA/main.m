Findex = 14;
dim = 100;
[low, up] = test_functions_range(Findex);
lb = low * ones(1, dim);
ub = up * ones(1, dim);
T = 1000;
populationSize = 30;
rng(1); % 'philox' is not supported in ga

options = optimoptions('ga', 'MaxGenerations', T, 'PopulationSize', populationSize);


% Define the objective function (test_functions)
objective = @(x) test_functions(x, Findex);

% Run Genetic Algorithm
[x, fval] = ga(objective, dim, [], [], [], [], lb, ub, [], options);

% Display results
bestsol = x;
bestfitness = fval;
disp(['Best solution: ', num2str(bestsol)]);
disp(['Best fitness: ', num2str(bestfitness)]);
