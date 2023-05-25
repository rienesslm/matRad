classdef NSGAOptimizer
    % matRad_DerivativeFreeOptimizer implements an derivative-free optimizer.
    %
    % This MATLAB code implements the NSGA (Non-dominated Sorting Genetic Algorithm)
    % for multi-objective optimization with an arbitrary number of objectives and constraints.
    % The `NSGA` class represents the main algorithm, while the `Individual` class represents
    % a solution with its objectives, constraints, rank, crowding distance, and domination count.
    %
    % The code follows a similar structure as the Java implementation.
    % It includes methods for initialization, crossover, dominance comparison,
    % cloning, and other required operations.
    % The `run` method executes the NSGA algorithm, and the `main` section demonstrates
    % an example usage.
    %
    % Please note that the code provided is a basic implementation of the NSGA in MATLAB
    % and may require further modifications and optimizations based on your specific
    % requirements and problem domain.
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2019 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties
        populationSize
        maxGenerations
        crossoverProbability
        mutationProbability
        numObjectives
        numConstraints
    end

    methods
        function obj = NSGAOptimizer(populationSize, maxGenerations, crossoverProbability, mutationProbability, numObjectives, numConstraints)
            obj.populationSize = populationSize;
            obj.maxGenerations = maxGenerations;
            obj.crossoverProbability = crossoverProbability;
            obj.mutationProbability = mutationProbability;
            obj.numObjectives = numObjectives;
            obj.numConstraints = numConstraints;
        end

        function run(obj)
            population = obj.initializePopulation();

            for generation = 1:obj.maxGenerations
                offspringPopulation = obj.createOffspringPopulation(population);
                combinedPopulation = obj.combinePopulations(population, offspringPopulation);
                fronts = obj.fastNonDominatedSort(combinedPopulation);
                selectedPopulation = obj.selection(fronts);

                % Perform necessary operations on the selected population
                % based on your specific problem domain

                % Example: Display the best individual in the current generation
                bestIndividual = selectedPopulation(1);
                disp(['Generation ', num2str(generation), ': Best individual objectives = ', num2str(bestIndividual.objectives)]);

                population = selectedPopulation;
            end
        end

        function population = initializePopulation(obj)
            population = Individual.empty(obj.populationSize, 0);
            for i = 1:obj.populationSize
                population(i) = Individual(obj.numObjectives, obj.numConstraints);
            end
        end

        function offspringPopulation = createOffspringPopulation(obj, population)
            offspringPopulation = Individual.empty();
            randomIndices = randperm(obj.populationSize);

            for i = 1:2:obj.populationSize
                parent1 = population(randomIndices(i));
                parent2 = population(randomIndices(i+1));

                if rand() < obj.crossoverProbability
                    [offspring1, offspring2] = parent1.crossover(parent2);
                    offspringPopulation = [offspringPopulation, offspring1, offspring2];
                else
                    offspringPopulation = [offspringPopulation, parent1.clone(), parent2.clone()];
                end
            end
        end

        function combinedPopulation = combinePopulations(obj, population, offspringPopulation)
            combinedPopulation = [population, offspringPopulation];
        end

        function fronts = fastNonDominatedSort(obj, population)
            % Perform fast non-dominated sorting on the population
            % Implement your non-dominated sorting algorithm here
            % Return the fronts as a cell array
            % Each cell represents a front and contains the indices of individuals in that front
            % Example: fronts = {front1Indices, front2Indices, ...}
            % Use the Individual.dominates method for dominance comparison
        end

        function selectedPopulation = selection(obj, fronts)
            % Perform selection based on crowding distance
            % Implement your selection method here
            % Return the selected population
            % Use the Individual.calculateCrowdingDistance method for crowding distance calculation
        end

        function obj = optimize(obj,w0,optiProb,dij,cst)
             funcs.objective         = @(x) optiProb.matRad_objectiveFunction(x,dij,cst);
             funcs.constraints       = @(x) optiProb.matRad_constraintFunctions(x,dij,cst); 
        end    
    end
end


% Example usage:
% populationSize = 100;
% maxGenerations = 100;
% crossoverProbability = 0.9;
% mutationProbability = 0.1;
% numObjectives = 2;
% numConstraints = 1;
%
% optimizer = NSGAOptimizer(populationSize, maxGenerations, crossoverProbability, mutationProbability, numObjectives, numConstraints);
% optimizer.run();
