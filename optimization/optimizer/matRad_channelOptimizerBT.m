classdef matRad_channelOptimizerBT
    properties
        populationSize=10; % populationSize size
        chromosomePairs=8; % number of pairs of chromosomes to be crossovered
        mutationNumber=6; % number chromosomes to be mutated
        numberGenerations=1; % Total number of generations
        Optimal_ChannelCombination = [];
        phantom;
        Min_fitness_value;
        populationMatrix;
        MAX_NUMBER_ACTIVE_CHANNELS = 20;
        MIN_NUMBER_ACTIVE_CHANNELS = 20;
        pln;
        dij;
        cst;
        stf;
        channelStability;

        % To plot the distribution of the fitness values in each
        % generation. (Show how values get closer to each other)
        fitnessValueMatrix;
    end

    methods
        function obj = matRad_channelOptimizerBT(phantom, populationSizeSize, generationSize, maxChannels, minChannels, activeNeedles)
            obj.phantom = phantom;
            obj.populationSize = populationSizeSize;
            obj.numberGenerations = generationSize;
            obj.MIN_NUMBER_ACTIVE_CHANNELS=minChannels;
            obj.MAX_NUMBER_ACTIVE_CHANNELS = maxChannels;
            obj.channelStability = "fixed";
            obj.pln.propStf.template.activeNeedles = activeNeedles;
            obj.populationMatrix = obj.initPopulation(obj.populationSize);
            [obj.pln, obj.dij, obj.cst, obj.stf] = obj.initPln();
            %P=obj.initpopulationSize(obj.populationSize);
            K=0;
            [x1,y1]=size(obj.populationMatrix);
            P1=0;
            obj.fitnessValueMatrix = zeros(obj.numberGenerations, x1+2*obj.chromosomePairs+obj.mutationNumber -1);
            for i=1:obj.numberGenerations
                Cr=obj.crossover(obj.populationMatrix, obj.chromosomePairs);
                Mu=obj.mutation(obj.populationMatrix, obj.mutationNumber);
                obj.populationMatrix(obj.populationSize + 1 : obj.populationSize + 2 * obj.chromosomePairs, :)=Cr;
                obj.populationMatrix(obj.populationSize + 2 * obj.chromosomePairs + 1 : obj.populationSize + 2 * obj.chromosomePairs + obj.mutationNumber, :)=Mu;
                objectiveValues=evaluation(obj, obj.populationMatrix);
                obj.fitnessValueMatrix(i, :) = objectiveValues;
                [obj.populationMatrix,S]=selection(obj, obj.populationMatrix, objectiveValues, obj.populationSize); %No clue why the function call needs 'obj' here...
                K(i,1)=sum(S)/obj.populationSize;
                K(i,2)=S(1); %best
            end
            obj.Min_fitness_value=min(K(:,2))
            obj.Optimal_ChannelCombination = obj.populationMatrix(1,:) % Best chromosome
            createBoxplotDiagram(obj.fitnessValueMatrix');
        end


        function populationMatrix = initPopulation(obj, pS)
            % Function creates a set of individuals. Each individual is
            % limited in the size of active
            % pS = populationSize size
            channelVectorSize=13*13;
%             baseChannels = [0 0 0 1 0 1 0 1 0 1 0 0 0;
%                             0 0 1 0 1 0 0 0 1 0 1 0 0;... % 6.5
%                             0 1 0 1 0 1 0 1 0 1 0 1 0;... % 6.0
%                             1 0 1 0 1 0 0 0 1 0 1 0 1;... % 5.5
%                             0 1 0 1 0 1 0 1 0 1 0 1 0;... % 5.0
%                             1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
%                             0 1 0 1 0 1 0 1 0 1 0 1 0;... % 4.0
%                             1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
%                             0 1 0 1 0 1 0 1 0 1 0 1 0;... % 3.0
%                             1 0 1 0 1 0 1 0 1 0 1 0 1;... % 2.5
%                             0 1 0 1 0 1 0 1 0 1 0 1 0;... % 2.0
%                             1 0 1 0 1 0 0 0 0 0 1 0 1;... % 1.5
%                             0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
%                            %A a B b C c D d E e F f G
            populationMatrix = zeros(pS + 1,channelVectorSize);
            A = [0	0	0	0	0	0	0	0	1	0	0	0	0
                0	0	0	0	0	0	0	0	0	0	0	0	1
                1	0	0	0	0	0	0	0	0	0	0	0	0
                0	0	0	0	0	0	0	0	0	0	0	0	0
                0	0	0	0	0	0	1	0	0	0	0	0	0
                1	0	0	0	0	0	0	0	0	0	0	0	0
                0	0	0	0	0	0	0	0	0	0	0	0	0
                0	0	0	0	0	0	0	0	0	0	0	1	0
                0	0	0	0	0	0	0	0	0	0	0	0	0
                0	0	0	0	0	0	0	0	0	0	0	0	0
                0	0	0	0	0	0	0	1	0	0	0	0	0
                0	0	0	1	0	0	0	0	0	0	0	0	0
                0	0	0	0	0	0	0	0	0	0	0	0	1];
            obj.pln.propStf.template.activeNeedles = A(:)';
            if isfield(obj.pln,'propStf')
                populationMatrix(1, :) = obj.pln.propStf.template.activeNeedles;
            else
                channelVector = zeros(1, channelVectorSize);
                numChannels = randi([obj.MIN_NUMBER_ACTIVE_CHANNELS, obj.MAX_NUMBER_ACTIVE_CHANNELS]);

                % % Generate random indices to set to 1
                indices = randperm(channelVectorSize, numChannels);
                
                % Set the selected indices to 1
                channelVector(indices) = 1;

                % Shuffle the vector randomly
                channelVector = channelVector(randperm(channelVectorSize));
                % Assign the vector to the matrix
                populationMatrix(1, :) = channelVector;
            end
%             populationMatrix(1, :) = reshape(baseChannels, 1, []);
            for k = 2 : pS
                channelVector = zeros(1, channelVectorSize);
                numChannels = randi([obj.MIN_NUMBER_ACTIVE_CHANNELS, obj.MAX_NUMBER_ACTIVE_CHANNELS]);

                % % Generate random indices to set to 1
                indices = randperm(channelVectorSize, numChannels);
                
                % Set the selected indices to 1
                channelVector(indices) = 1;

                % Shuffle the vector randomly
                channelVector = channelVector(randperm(channelVectorSize));

                % Assign the vector to the matrix
                populationMatrix(k, :) = channelVector;

                % channelMatrix=round(rand(13,13));
                % channelVector=reshape(channelMatrix.',1,[]);
                % populationMatrix(k,:) = channelVector;
            end
            
            size(populationMatrix)
        end

        function [pln, dij, cst, stf] = initPln(obj)
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
            cst{5,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(15,20));
            cst{6,5}.Priority = 1;

            % In this example, the lymph nodes will not be part of the treatment:
            cst{7,6}    =  [];
            cst{7,3}    =  'OAR';

            %A PTV is not needed, but we will use it to control the dose fall off
            cst{6,3}    =  'OAR';
            cst{6,6}{1} =  struct(DoseObjectives.matRad_SquaredOverdosing(25,12));
            cst{6,5}.Priority = 2;

            % Rectum Objective
            cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(20,10));

            % Bladder Objective
            cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(40,14));

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

            pln.propDoseCalc.doseGrid.resolution.x = 10; % [mm]
            pln.propDoseCalc.doseGrid.resolution.y = 10; % [mm]
            pln.propDoseCalc.doseGrid.resolution.z = 10; % [mm]

            pln.propDoseCalc.DistanceCutoff    = 70; %[mm] sets the maximum distance
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

            dij = matRad_calcBrachyDose(ct,stf,pln,cst);
        end

        function Y=crossover(obj, populationMatrix, n)
            % populationMatrix = populationMatrix
            % n = number of pairs of chromosomes to be crossovered
            [x1,y1]=size(populationMatrix);
            Z=zeros(2*n,y1);
            if obj.channelStability == "fixed"
                i = 1;
                while i <= n
                    r1=randi(x1,1,2);
                    while r1(1)==r1(2)
                        r1=randi(x1,1,2);
                    end
                    A1=populationMatrix(r1(1),:); % parent 1
                    A2=populationMatrix(r1(2),:); % parent 2
                    r2=1+randi(y1-1); % random cutting point
                    B1=A1(1,r2:y1);
                    A1(1,r2:y1)=A2(1,r2:y1);
                    A2(1,r2:y1)=B1;
                    if nnz(A1) == obj.MAX_NUMBER_ACTIVE_CHANNELS && nnz(A2) == obj.MAX_NUMBER_ACTIVE_CHANNELS
                        Z(2*i-1,:)=A1; % offspring 1
                        Z(2*i,:)=A2; % offspring 2
                        i = i +1;
                    end
                end
            else
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
            end
            Y=Z;
        end

        function Y=mutation(obj, populationMatrix, n)
            % populationMatrix = populationMatrix
            % n = chromosomes to be mutated
            [x1,y1]=size(populationMatrix);
            Z=zeros(n,y1);
            if obj.channelStability == "fixed"
                i = 1;
                while i <= n
                    j = 0;
                    while j == 0
                        randomIndividualIndex=randi(x1);
                        randomParent=populationMatrix(randomIndividualIndex,:); % random parent
                        randomMutationPoint1=randi(y1);
                        randomMutationPoint2=randi(y1);
                        if randomParent(1,randomMutationPoint1) ~= randomParent(1,randomMutationPoint2)
                            j = 1;
                        end
                    end
                    if randomParent(1,randomMutationPoint1)== 1 && randomParent(1,randomMutationPoint2)== 0
                        randomParent(1,randomMutationPoint1) = 0; % flip the first bit
                        randomParent(1, randomMutationPoint2) = 1; % flip the second bit
                    elseif randomParent(1,randomMutationPoint1)== 0 && randomParent(1,randomMutationPoint2)== 1
                        randomParent(1,randomMutationPoint1) = 1; % flip the first bit
                        randomParent(1, randomMutationPoint2) = 0; % flip the second bit
                    end
                    if nnz(randomParent) == obj.MAX_NUMBER_ACTIVE_CHANNELS
                        Z(i,:)=randomParent;
                        i = i + 1;
                    end
                end
            else
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
            end
            Y=Z;
        end

        function [listOfFunctionValues]=evaluation(obj, populationMatrix)
            [x1, y1]=size(populationMatrix)
            listOfFunctionValues=zeros(1,x1);
            for i = 1:x1
                matRad_rc;
                load(obj.phantom);
                populationVector = populationMatrix(i,:);
                obj.pln.propStf.template.activeNeedles = reshape(populationVector, 13, 13);
                pln = obj.pln;
                %populationVector = populationMatrix(i,:);
                %pln.propStf.template.activeNeedles = reshape(populationVector, 13, 13);
                dij = obj.dij;
                cst = obj.cst;
                stf = obj.stf;
                dwellPositions = obj.getDwellPositions(cst, pln);
                dij.dwellPositions = dwellPositions;
                dij.state = "Genetic Optimze";
                %populationVector = populationMatrix(i,:);
                %pln.propStf.template.activeNeedles = reshape(populationVector, 13, 13);
                H = matRad_fluenceOptimization(dij, cst, pln);
                listOfFunctionValues(1,i) = H.info.objective;
            end
            listOfFunctionValues
        end

        function [newPopulationMatrix, fitnessValues] = selection(obj, populationMatrix, F, pS)
            % populationMatrix - populationMatrix, F -fitness value, pS - population size
            [x,y]=size(populationMatrix);
            Y1=zeros(pS,y);
            F = F +10; % adding 10 to ensure no chromosome has negative fitness
            % elite selection
            e=3;
            for i = 1:e
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
            while d2 <= pS-1
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
            newPopulationMatrix=Y1;
            fitnessValues=Fn - 10; % substract 10 to return the original fitness
        end

        function decimalNumber = matRad_bit2de(obj, binaryString)
            % Convert binary strings to decimal numbers
            decimalNumber = polyval(binaryString, 2);
        end

        function selectedDwellPositions = getDwellPositions(obj, cst, pln)
            load(obj.phantom);
            %% config
            matRad_cfg = MatRad_Config.instance();
            addpath(fullfile( matRad_cfg.matRadRoot));
            matRad_cfg.dispInfo('matRad: Generating stf struct... ');

            if ~isfield(pln,'propStf')
                matRad_cfg.dispError('no applicator information in pln struct');
            end



            %% generate image coordinates

            % find all target voxels from cst cell array
            V = [];
            for i=1:size(cst,1)
                if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
                    V = [V;vertcat(cst{i,4}{:})];
                end
            end

            % Remove double voxels
            V = unique(V);
            % generate voi cube for targets
            voiTarget    = zeros(ct.cubeDim);
            voiTarget(V) = 1;

            % add margin
            addmarginBool = matRad_cfg.propStf.defaultAddMargin;
            if isfield(pln,'propStf') && isfield(pln.propStf,'addMargin')
                addmarginBool = pln.propStf.addMargin;
            end

            if addmarginBool
                voiTarget = matRad_addMargin(voiTarget,cst,ct.resolution,ct.resolution,true);
                V   = find(voiTarget>0);
            end

            % throw error message if no target is found
            if isempty(V)
                matRad_cfg.dispError('Could not find target.');
            end

            % Convert linear indices to 3D voxel coordinates
            [coordsY_vox, coordsX_vox, coordsZ_vox] = ind2sub(ct.cubeDim,V);

            %translate to geometric coordinates and save in stf

            targetVolume.Xvox = ct.x(coordsX_vox); % angabe in mm
            targetVolume.Yvox = ct.y(coordsY_vox);
            targetVolume.Zvox = ct.z(coordsZ_vox);


            % % calculate rED or rSP from HU
            % ct = matRad_calcWaterEqD(ct, pln);

            % take only voxels inside patient
            V = [cst{:,4}];
            V = unique(vertcat(V{:}));

            % ignore densities outside of contours
            eraseCtDensMask = ones(prod(ct.cubeDim),1);
            eraseCtDensMask(V) = 0;
            for i = 1:ct.numOfCtScen
                ct.cube{i}(eraseCtDensMask == 1) = 0;
            end

            %% save VOI coordinates
            %translate to geometric coordinates and save in stf
            targetVolume.Xvox = ct.x(coordsX_vox);
            targetVolume.Yvox = ct.y(coordsY_vox);
            targetVolume.Zvox = ct.z(coordsZ_vox);
            %         stf.targetVolume.Xvox = ct.x(coordsX_vox); % given in mm
            %         stf.targetVolume.Yvox = ct.y(coordsY_vox);
            %         stf.targetVolume.Zvox = ct.z(coordsZ_vox);
            %% meta info from pln
            radiationMode = pln.radiationMode;
            numOfSeedsPerNeedle = pln.propStf.needle.seedsNo;
            numOfNeedles = nnz(pln.propStf.template.activeNeedles);
            totalNumOfBixels = numOfSeedsPerNeedle*numOfNeedles; % means total number of seeds

            %% generate 2D template points
            % the template origin is set at its center. In the image coordinate system,
            % the center will be positioned at the bottom of the volume of interest.

            [row,col] = find(pln.propStf.template.activeNeedles);
            templX = col*pln.propStf.bixelWidth + pln.propStf.templateRoot(1) - (13+1)/2*pln.propStf.bixelWidth;
            templY = row*pln.propStf.bixelWidth + pln.propStf.templateRoot(2) - (13+1)/2*pln.propStf.bixelWidth;
            templZ = ones(size(col))                 + pln.propStf.templateRoot(3);


            template = [templX';templY';templZ'];




            %% generate seed positions
            % seed positions can be generated from neeldes, template and oriantation
            % needles are assumed to go trough the template vertically

            % needle position
            d = pln.propStf.needle.seedDistance;
            seedsNo = pln.propStf.needle.seedsNo;
            needleDist(1,1,:) = d.*[0:seedsNo-1]'; % 1x1xN Array with seed positions on needle
            needleDir = needleDist.*[0;0;1];
            seedPos_coord_need_seed = needleDir + template;
            seedPos_need_seed_coord = shiftdim(seedPos_coord_need_seed,1);
            % the output array has the dimentions (needleNo,seedNo,coordinates)
            X = seedPos_need_seed_coord(:,:,1);
            Y = seedPos_need_seed_coord(:,:,2);
            Z = seedPos_need_seed_coord(:,:,3);

            seedPoints.x = reshape(X,1,[]);
            seedPoints.y = reshape(Y,1,[]);
            seedPoints.z = reshape(Z,1,[]);
            selectedDwellPositions = zeros(obj.stf.totalNumOfBixels, 1);
            
            % ind = find(obj.stf.seedPoints.x == seedPoints.x && obj.stf.seedPoints.x == seedPoints.y && obj.stf.seedPoints.z == seedPoints.z);

            for f = 1:obj.stf.totalNumOfBixels
                for j = 1:totalNumOfBixels
                    if isequal(obj.stf.seedPoints.x(f), seedPoints.x(j)) && isequal(obj.stf.seedPoints.y(f), seedPoints.y(j)) && isequal(obj.stf.seedPoints.z(f), seedPoints.z(j))
                        selectedDwellPositions(f) = 1;
                    end
                end
            end
            % ind = find(obj.stf.seedPoints.x == seedPoints.x && obj.stf.seedPoints.x == seedPoints.y && obj.stf.seedPoints.z == seedPoints.z);
            selectedDwellPositions;
        end
    end
end
