classdef matRad_channelOptimizerBT
    properties     
        populationSize=10; % populationSize size
        chromosomePairs=2; % number of pairs of chromosomes to be crossovered
        mutationNumber=2; % number chromosomes to be mutated
        numberGenerations=5; % Total number of generations
        Optimal_ChannelCombination = [];
        phantom;
        Min_fitness_value;
        populationMatrix;
        pln;
        dij;
        cst;
    end
    
    methods
        function obj = matRad_channelOptimizerBT(phantom, populationSizeSize, generationSize)
            obj.phantom = phantom;
            obj.populationSize = populationSizeSize;
            obj.numberGenerations = generationSize;
            obj.populationMatrix = obj.initpopulationSize(obj.populationSize);
            obj.pln = obj.initPln();
            %P=obj.initpopulationSize(obj.populationSize);
            K=0;
            [x1,y1]=size(obj.populationMatrix);
            P1=0;
            for i=1:obj.numberGenerations
                Cr=obj.crossover(obj.populationMatrix, obj.chromosomePairs);
                Mu=obj.mutation(obj.populationMatrix, obj.mutationNumber);
                obj.populationMatrix(obj.populationSize + 1 : obj.populationSize + 2 * obj.chromosomePairs, :)=Cr;
                obj.populationMatrix(obj.populationSize + 2 * obj.chromosomePairs + 1 : obj.populationSize + 2 * obj.chromosomePairs + obj.mutationNumber, :)=Mu;
                E=evaluation(obj, obj.populationMatrix)
                [obj.populationMatrix,S]=selection(obj, obj.populationMatrix, E, obj.populationSize); %No clue why the function call needs 'obj' here...
                K(i,1)=sum(S)/obj.populationSize;
                K(i,2)=S(1); %best
            end
            obj.Min_fitness_value=min(K(:,2))
            obj.Optimal_ChannelCombination = obj.populationMatrix(1,:); % Best chromosome
        end
        

        function populationMatrix = initpopulationSize(obj, pS)
            % pS = populationSize size
            % Y = cell(pS, 1) ;
            populationMatrix = zeros(pS,13*13);
            for k = 1 : pS
                channelMatrix=round(rand(13,13));
                channelVector=reshape(channelMatrix.',1,[]);
                populationMatrix(k,:) = channelVector;
            end
            size(populationMatrix)
        end

        function pln = initPln(obj)
            matRad_rc;
                load(obj.phantom);

                %% I - update/set dose objectives for brachytherapy
                % The sixth column represents dose objectives and constraints for
                % optimization: First, the objective function for the individual structure
                % is chosen, its parameters denote doses that should not be tranceeded
                % towards higher or lower doses (SquaredOverdose, SquaredUnderdose) or
                % doses that are particularly aimed for (SquaredUnderDose).

                disp(cst{6,6}{1});

                % Following frequently prescribed planning doses of 15 Gy
                % (https://pubmed.ncbi.nlm.nih.gov/22559663/) objectives can be updated to:

                % the planning target was changed to the clinical segmentation of the
                % prostate bed.
                cst{5,3}    = 'TARGET';
                cst{5,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(100,15));
                cst{5,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(100,17.5));
                cst{6,5}.Priority = 1;

                % In this example, the lymph nodes will not be part of the treatment:
                cst{7,6}    =  [];
                cst{7,3}    =  'OAR';

                %A PTV is not needed, but we will use it to control the dose fall off
                cst{6,3}    =  'OAR';
                cst{6,6}{1} =  struct(DoseObjectives.matRad_SquaredOverdosing(100,12));
                cst{6,5}.Priority = 2;

                % Rectum Objective
                cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,7.5));

                % Bladder Objective
                cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,7.5));

                % Body NT objective
                cst{9,6}{1} = struct(DoseObjectives.matRad_MeanDose(1));


                %% II.1 Treatment Plan
                % The next step is to define your treatment plan labeled as 'pln'. This
                % matlab structure requires input from the treatment planner and defines
                % the most important cornerstones of your treatment plan.

                % First of all, we need to define what kind of radiation modality we would
                % like to use. Possible values are photons, protons or carbon
                % (external beam) or brachy as an invasive tratment option.
                % In this case we want to use brachytherapy. Then, we need to define a
                % treatment machine to correctly load the corresponding base data.
                % matRad includes example base data for HDR and LDR brachytherapy.
                % Here we will use HDR. By this means matRad will look for 'brachy_HDR.mat'
                % in our root directory and will use the data provided in there for
                % dose calculation.

                pln.radiationMode   = 'brachy';
                pln.machine         = 'HDR';    % 'LDR' or 'HDR' for brachy


                %% II.1 - needle and template geometry
                % Now we have to set some parameters for the template and the needles.
                % Let's start with the needles: Seed distance is the distance between
                % two neighbouring seeds or holding points on one needle or catheter. The
                % seeds No denotes how many seeds/holding points there are per needle.

                pln.propStf.needle.seedDistance      = 10; % [mm]
                pln.propStf.needle.seedsNo           = 6;

                %% II.1 - template position
                % The implantation is normally done through a 13 x 13 template from the
                % patients inferior, which is the negative z axis here.
                % The direction of the needles is defined by template normal.
                % Neighbour distances are called by bixelWidth, because this field is also
                % used for external beam therapy.
                % The needles will be positioned right under the target volume pointing up.

                pln.propStf.template.normal      = [0,0,1];
                pln.propStf.bixelWidth   = 5; % [mm] template grid distance
                pln.propStf.templateRoot = matRad_getTemplateRoot(ct,cst); % mass center of
                % target in x and y and bottom in z

                % Here, we define active needles as 1 and inactive needles
                % as 0. This is the x-y plane and needles point in z direction.
                % A checkerboard pattern is frequantly used. The whole geometry will become
                % clearer when it is displayed in 3D view in the next section.

                % Initial channel distribution for precalculation of dij.
                pln.propStf.template.activeNeedles =   [1 1 1 1 1 1 1 1 1 1 1 1 1;... % 7.0
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 6.5
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 6.0
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 5.5
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 5.0
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 4.5
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 4.0
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 4.5
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 3.0
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 2.5
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 2.0
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1;... % 1.5
                                                        1 1 1 1 1 1 1 1 1 1 1 1 1];   % 1.0
                                                       %A a B b C c D d E e F f G
                
                pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0); %  target center
                
                %% II.1 - dose calculation options
                % for dose calculation we use eather the 2D or the 1D formalism proposed by
                % TG 43. Also, set resolution of dose calculation and optimization.
                % If your system gets stuck with the resolution, you can lower it to 10 or
                % 20, just to get an initial result. Otherwise, reduce the number of
                % needles.
                % Calculation time will be reduced by one tenth when we define a dose
                % cutoff distance.
                pln.propDoseCalc.TG43approximation = '2D'; %'1D' or '2D'

                pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
                pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
                pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

                pln.propDoseCalc.DistanceCutoff    = 130; %[mm] sets the maximum distance
                %to which dose is calculated.

                % the standard interior point optimizer IPOPT can be used
                pln.propOpt.optimizer       = 'IPOPT';

                %% II.1 - book keeping
                % Some field names have to be kept although they don't have a direct
                % relevance for brachy therapy.
                pln.propOpt.bioOptimization = 'none';
                pln.propOpt.runDAO          = false;
                pln.propOpt.runSequencing   = false;
                pln.propStf.gantryAngles    = [];
                pln.propStf.couchAngles     = [];
                pln.propStf.numOfBeams      = 0;
                pln.numOfFractions          = 1;

                 %% II.1 - view plan
                % Et voila! Our treatment plan structure is ready. Lets have a look:
                %disp(pln);


                %% II.2 Steering Seed Positions From STF
                % The steering file struct contains all needls/catheter geometry with the
                % target volume, number of needles, seeds and the positions of all needles
                % The one in the end enables visualization.

                stf = matRad_generateStf(ct,cst,pln,1);

                %% II.2 - view stf
                % The 3D view is interesting, but we also want to know how the stf struct
                % looks like.

                %disp(stf)

                %% II.3 - Dose Calculation
                % Let's generate dosimetric information by pre-computing a dose influence
                % matrix for seed/holding point intensities. Having dose influences
                % available allows subsequent inverse optimization.
                % Don't get inpatient, this can take a few seconds...

                obj.dij = matRad_calcBrachyDose(ct,stf,pln,cst);
                obj.cst = cst;
                
        end

        function Y=crossover(obj, populationMatrix, n)
            % populationMatrix = populationMatrix
            % n = number of pairs of chromosomes to be crossovered
            [x1,y1]=size(populationMatrix);
            Z=zeros(2*n,y1);
            for i = 1:n
                r1=randi(x1,1,2);
                while r1(1)==r1(2)
                    r1=randi(x1,1,2);
                end
                A1=populationMatrix(r1(1),:); % parent 1
                A2=populationMatrix(r1(2),:); % parent 2
                r2=1+randi(y1-1); % random cutting point
                B1=A1(1,r2:y1);
                A1(1,r2:y1)=A2(1,r2:y1);
                %A1(1,r2:y1)=A2(1,r2:9);
                A2(1,r2:y1)=B1;
                %A2(1,r2:9)=B1;
                Z(2*i-1,:)=A1; % offspring 1
                Z(2*i,:)=A2; % offspring 2
            end
            Y=Z;
        end

        function Y=mutation(obj, populationMatrix, n)
            % populationMatrix = populationMatrix
            % n = chromosomes to be mutated
            [x1,y1]=size(populationMatrix);
            Z=zeros(n,y1);
            for i = 1:n
                r1=randi(x1);
                A1=populationMatrix(r1,:); % random parent
                r2=randi(y1);
                if A1(1,r2)== 1
                    A1(1,r2) = 0; % flip the bit
                else
                    A1(1,r2) = 1;
                end
                Z(i,:)=A1;
            end
            Y=Z;
        end

        function [listOfFunctionValues]=evaluation(obj, populationMatrix)
            [x1, y1]=size(populationMatrix)
            listOfFunctionValues=zeros(1,x1);
            for i = 1:x1                
                populationVector = populationMatrix(i,:);
                obj.pln.propStf.template.activeNeedles = reshape(populationVector, 13, 13);
                %% III Inverse Optimization for brachy therapy
                % The goal of the fluence optimization is to find a set of holding point
                % times which yield the best possible dose distribution according to
                % the clinical objectives and constraints underlying the radiation
                % treatment. Once the optimization has finished, trigger to
                % visualize the optimized dose cubes.                
                H = matRad_fluenceOptimization(obj.dij, obj.cst, obj.pln);
                listOfFunctionValues(1,i) = H.info.objective;
            end
            listOfFunctionValues
        end

        function [YY1, YY2] = selection(obj, populationMatrix, F, pS)
            % populationMatrix - populationMatrix, F -fitness value, pS - population size
            [x,y]=size(populationMatrix);
            Y1=zeros(pS,y);
            %F = F +10; % adding 10 to ensure no chromosome has negative fitness
            % elite selection
            e=3;
            for i = 1:3
                obj.Min_fitness_value = min(F);
                find(F==min(F))
                [r1,c1]=find(F==min(F))
                populationMatrix(min(c1),:)
                Y1(i,:)=populationMatrix(min(c1),:);
                populationMatrix(min(c1),:)=[];
                Fn(i)=F(min(c1));
                F(:,min(c1))=[];
            end
            D=F/sum(F); % Determine selection probability
            E=cumsum(D); % Determine cumulative probability
            N=rand(1); % Generate a vector constaining normalised random numbers
            d1=1;
            d2=e;
            while d2 <= pS-e
                if N <= E(d1)
                    Y1(d2+1,:)=populationMatrix(d1,:);
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

