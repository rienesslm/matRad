classdef beamAngleOptimizer
    properties     
        population=10; % population size
        chromosomePairs=2; % number of pairs of chromosomes to be crossovered
        mutationNumber=2; % number chromosomes to be mutated
        numberGenerations=10; % Total number of generations
    end
    
    methods
        function obj = beamAngleOptimizer()
            P=obj.initPopulation(obj.population)
            K=0;
            [x1,y1]=size(P);
            P1=0;
            for i=1:obj.numberGenerations
                Cr=obj.crossover(P, obj.chromosomePairs);
                Mu=obj.mutation(P, obj.mutationNumber);
                P(obj.population + 1 : obj.population + 2 * obj.chromosomePairs, :)=Cr;
                P(obj.population + 2 * obj.chromosomePairs + 1 : obj.population + 2 * obj.chromosomePairs + obj.mutationNumber, :)=Mu;
                E=evaluation(obj, P)
                [P,S]=selection(obj, P, E, obj.population); %No clue why the function call needs 'obj' here...
                K(i,1)=sum(S)/obj.population;
                K(i,2)=S(1); %best
            end
            Min_fitness_value=min(K(:,2))
            P2 = P(1,:); % Best chromosome
            % convert binary to real number
            
            A=mod(obj.matRad_bit2de(P2(1, 1:y1/2)), 360);
            B=mod(obj.matRad_bit2de(P2(1, y1/2+1:y1)), 360);
            Optimal_Angels=[A,B]
        end
        

        function Y = initPopulation(obj, n)
            % n = population size

            % It is noted that the number of bits to represent the variables
            % in binary numbers depends on the required accuracy (the number
            % of digits after comma)

            % In this example, I want the solution precision with 5 places after the
            % decimal point, and with the upper and lower bounds of the variables are 3
            % and -3, so, for each variable, we need 20 bits.
            % General formula: 2^(m-1) < (upper bound - lower bound)*10^p < 2^m -1
            % In this case: p = 5 and m = 20.

            % We have 2 variables (x and y), so we need 40 bits in
            % total for binary encoding
            %Y=round(rand(n,40));
            Y=round(rand(n,18));

        end

        function Y=crossover(obj, P, n)
            % P = population
            % n = number of pairs of chromosomes to be crossovered
            [x1,y1]=size(P);
            Z=zeros(2*n,y1);
            for i = 1:n
                r1=randi(x1,1,2);
                while r1(1)==r1(2)
                    r1=randi(x1,1,2);
                end
                A1=P(r1(1),:); % parent 1
                A2=P(r1(2),:); % parent 2
                r2=1+randi(y1-1); % random cutting point
                B1=A1(1,r2:y1);
                A1(1,r2:y1)=A2(1,r2:18);
                A2(1,r2:18)=B1;
                Z(2*i-1,:)=A1; % offspring 1
                Z(2*i,:)=A2; % offspring 2
            end
            Y=Z;
        end

        function Y=mutation(obj, P, n)
            % P = population
            % n = chromosomes to be mutated
            [x1,y1]=size(P);
            Z=zeros(n,y1);
            for i = 1:n
                r1=randi(x1);
                A1=P(r1,:); % random parent
                r2=randi(y1);
                if A1(1,r2)== 1
                    A1(1,r2) = 0; % flick the bit
                else
                    A1(1,r2) = 1;
                end
                Z(i,:)=A1;
            end
            Y=Z;
        end

        function [listOfFunctionValues]=evaluation(obj, P)
            [x1, y1]=size(P)
            listOfFunctionValues=zeros(1,x1);
            for i = 1:x1
                % (vi) how to compare the two results

                %% Patient Data Import
                % Let's begin with a clear Matlab environment and import the prostate
                % patient into your workspace

                matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

                load('TG119.mat');

                %% Treatment Plan
                % The next step is to define your treatment plan labeled as 'pln'. This
                % structure requires input from the treatment planner and defines the most
                % important cornerstones of your treatment plan.

                %%
                % First of all, we need to define what kind of radiation modality we would
                % like to use. Possible values are photons, protons or carbon. In this
                % example we would like to use protons for treatment planning. Next, we
                % need to define a treatment machine to correctly load the corresponding
                % base data. matRad features generic base data in the file
                % 'proton_Generic.mat'; consequently the machine has to be set accordingly
                pln.radiationMode = 'protons';
                pln.machine       = 'Generic';

                %%
                % Define the flavor of biological optimization for treatment planning along
                % with the quantity that should be used for optimization. Possible values
                % are (none: physical optimization; const_RBExD: constant RBE of 1.1;
                % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of
                % RBE-weighted dose. As we use protons, we follow here the clinical
                % standard and use a constant relative biological effectiveness of 1.1.
                % Therefore we set bioOptimization to const_RBExD
                pln.propOpt.bioOptimization = 'const_RBExD';

                %%
                % for particles it is possible to also calculate the LET disutribution
                % alongside the physical dose. Therefore you need to activate the
                % corresponding option during dose calculcation
                pln.propDoseCalc.calcLET = 1;

                %%
                % Now we have to set the remaining plan parameters.
                pln.numOfFractions        = 30;
                pln.propStf.gantryAngles  = [90 270];
                pln.propStf.couchAngles   = [0 0];
                pln.propStf.bixelWidth    = 3;
                pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
                pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
                pln.propOpt.runDAO        = 0;
                pln.propOpt.runSequencing = 0;

                % dose calculation settings
                pln.propDoseCalc.doseGrid.resolution.x = 10; % [mm]
                pln.propDoseCalc.doseGrid.resolution.y = 10; % [mm]
                pln.propDoseCalc.doseGrid.resolution.z = 10; % [mm]
                P
                IndividumA = P(i,1:y1/2)
                IndividumB = P(i,y1/2+1:y1)
                %IndividumA = P(i:x1/2,:)
                %IndividumB = P(x1/2+1:x1,:)
                %valueIndividumA = obj.matRad_bit2de(P(i,1:y1/2))
                valueIndividumA = obj.matRad_bit2de(IndividumA)
                %valueIndividumB = obj.matRad_bit2de(P(i,y1/2+1:y1))
                valueIndividumB = obj.matRad_bit2de(IndividumB)
                %alpha=mod(obj.matRad_bit2de(P(i,1:y1/2)), 360)
                alpha=mod(valueIndividumA, 360)
                %beta=mod(obj.matRad_bit2de(P(i,y1/2+1:y1)), 360)
                beta=mod(valueIndividumB, 360)
                pln.propStf.gantryAngles = [alpha beta];
                %% Generate Beam Geometry STF
                stf = matRad_generateStf(ct,cst,pln);

                %% Dose Calculation
                % Lets generate dosimetric information by pre-computing dose influence
                % matrices for unit beamlet intensities. Having dose influences available
                % allows for subsequent inverse optimization.
                dij = matRad_calcParticleDose(ct,stf,pln,cst);

                %% Inverse Optimization for IMPT
                % The goal of the fluence optimization is to find a set of bixel/spot
                % weights which yield the best possible dose distribution according to the
                % clinical objectives and constraints underlying the radiation treatment
                H = matRad_fluenceOptimization(dij,cst,pln);
                listOfFunctionValues(1,i) = H.info.objective;
            end
            listOfFunctionValues
        end

        function [YY1, YY2] = selection(obj, P, F, p)
            % P - population, F -fitness value, p - population size
            [x,y]=size(P);
            Y1=zeros(p,y);
            %F = F +10; % adding 10 to ensure no chromosome has negative fitness
            % elite selection
            e=3;
            for i = 1:3
                find(F==min(F))
                [r1,c1]=find(F==min(F))
                P(min(c1),:)
                Y1(i,:)=P(min(c1),:);
                P(min(c1),:)=[];
                Fn(i)=F(min(c1));
                F(:,min(c1))=[];
            end
            D=F/sum(F); % Determine selection probability
            E=cumsum(D); % Determine cumulative probability
            N=rand(1); % Generate a vector constaining normalised random numbers
            d1=1;
            d2=e;
            while d2 <= p-e
                if N <= E(d1)
                    Y1(d2+1,:)=P(d1,:);
                    Fn(d2+1)=F(d1);
                    N=rand(1);
                    d2=d2+1;
                    d1=1;
                else
                    d1 = d1 + 1;
                end
            end
            YY1=Y1;
            YY2=Fn; % substract 10 to return the original fitness
        end
        
        function decimalNumber = matRad_bit2de(obj, binaryString)
            % Convert binary strings to decimal numbers
            decimalNumber = polyval(binaryString, 2);
        end
    end
end

