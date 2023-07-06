classdef matRad_channelOptimizerBT
    properties     
        population=10; % population size
        chromosomePairs=2; % number of pairs of chromosomes to be crossovered
        mutationNumber=2; % number chromosomes to be mutated
        numberGenerations=5; % Total number of generations
        Optimal_ChannelCombination = [];
        phantom;
        Min_fitness_value;
    end
    
    methods
        function obj = matRad_channelOptimizerBT(phantom, populationSize, generationSize)
            obj.phantom = phantom;
            obj.population = populationSize;
            obj.numberGenerations = generationSize;
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
            obj.Min_fitness_value=min(K(:,2))
            P2 = P(1,:); % Best chromosome
            % convert binary to real number

            for i=1:obj.beamNumber
                angle=mod(obj.matRad_bit2de(P2(1, 1+(9*(i-1)):y1*i/obj.beamNumber)), 360);
                obj.Optimal_Angels(i)=angle;
            end
            %A=mod(obj.matRad_bit2de(P2(1, 1:y1/2)), 360);
            %B=mod(obj.matRad_bit2de(P2(1, y1/2+1:y1)), 360);
            %A=mod(obj.matRad_bit2de(P2(1, 1:y1)), 360);
            %Optimal_Angels = [A B]
            obj.Optimal_Angels
        end
        

        function Y = initPopulation(obj, pS)
            % pS = population size
            Y = cell(pS, 1) ;
            for k = 1 : Ps
                Y{k} = round(rand(13,13)) ;
            end

            %Y=round(rand(13,13));

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
                A1(1,r2:y1)=A2(1,r2:13);
                %A1(1,r2:y1)=A2(1,r2:9);
                A2(1,r2:13)=B1;
                %A2(1,r2:9)=B1;
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
                    A1(1,r2) = 0; % flip the bit
                else
                    A1(1,r2) = 1;
                end
                Z(i,:)=A1;
            end
            Y=Z;
        end

        function [listOfFunctionValues]=evaluation(obj, P)
            x1 = size(P)
            listOfFunctionValues=zeros(1,x1);
            for i = 1:x1
                %% List of contents
                % In this example we will show
                % (i) how to load patient data into matRad
                % (ii) how to setup an HDR brachy dose calculation and
                % (iii) how to inversely optimize holding position intensties
                % (iv) how to visually and quantitatively evaluate the result
                % (v) how to verify that functions do the right thing

                %% I Patient Data Import
                % Let's begin with a clear Matlab environment. Then, import the TG119
                % phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
                % structure defining the CT images and the structure set. Make sure the
                % matRad root directory with all its subdirectories is added to the Matlab
                % search path.

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

%                 pln.propStf.template.activeNeedles =   [0 0 0 1 0 1 0 1 0 1 0 0 0;... % 7.0
%                                                         0 0 1 0 1 0 0 0 1 0 1 0 0;... % 6.5
%                                                         0 1 0 1 0 1 0 1 0 1 0 1 0;... % 6.0
%                                                         1 0 1 0 1 0 0 0 1 0 1 0 1;... % 5.5
%                                                         0 1 0 1 0 1 0 1 0 1 0 1 0;... % 5.0
%                                                         1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
%                                                         0 1 0 1 0 1 0 1 0 1 0 1 0;... % 4.0
%                                                         1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
%                                                         0 1 0 1 0 1 0 1 0 1 0 1 0;... % 3.0
%                                                         1 0 1 0 1 0 1 0 1 0 1 0 1;... % 2.5
%                                                         0 1 0 1 0 1 0 1 0 1 0 1 0;... % 2.0
%                                                         1 0 1 0 1 0 0 0 0 0 1 0 1;... % 1.5
%                                                         0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
%                                                        %A a B b C c D d E e F f G
                
                pln.propStf.template.activeNeedles = P(i);

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
                disp(pln);


                %% II.2 Steering Seed Positions From STF
                % The steering file struct contains all needls/catheter geometry with the
                % target volume, number of needles, seeds and the positions of all needles
                % The one in the end enables visualization.

                stf = matRad_generateStf(ct,cst,pln,1);

                %% II.2 - view stf
                % The 3D view is interesting, but we also want to know how the stf struct
                % looks like.

                disp(stf)

                %% II.3 - Dose Calculation
                % Let's generate dosimetric information by pre-computing a dose influence
                % matrix for seed/holding point intensities. Having dose influences
                % available allows subsequent inverse optimization.
                % Don't get inpatient, this can take a few seconds...

                dij = matRad_calcBrachyDose(ct,stf,pln,cst);

                %% III Inverse Optimization for brachy therapy
                % The goal of the fluence optimization is to find a set of holding point
                % times which yield the best possible dose distribution according to
                % the clinical objectives and constraints underlying the radiation
                % treatment. Once the optimization has finished, trigger to
                % visualize the optimized dose cubes.

                resultGUI = matRad_fluenceOptimization(dij,cst,pln);
                % (vi) how to compare the two results

                %% Patient Data Import
                % Let's begin with a clear Matlab environment and import the prostate
                % patient into your workspace

                matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
                
                load(obj.phantom);
                %load('LIVER.mat');

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
                pln.propStf.gantryAngles  = zeros(1,obj.beamNumber);
                %pln.propStf.gantryAngles  = [90];
                pln.propStf.couchAngles   = zeros(1,obj.beamNumber);
                %pln.propStf.couchAngles   = [0];
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
                %IndividumA = P(i,1:y1/2)
                %(1, 1+(9*(i-1)):y1*i/obj.beamNumber)
                for j=1:obj.beamNumber
                    angle=mod(obj.matRad_bit2de(P(i, 1+(9*(j-1)):y1*j/obj.beamNumber)), 360);
                    pln.propStf.gantryAngles(j)=angle;
                end
                %IndividumA = P(i,1:y1)
                %IndividumB = P(i,y1/2+1:y1)
                %IndividumA = P(i:x1/2,:)
                %IndividumB = P(x1/2+1:x1,:)
                %valueIndividumA = obj.matRad_bit2de(P(i,1:y1/2))
                %valueIndividumA = obj.matRad_bit2de(IndividumA)
                %valueIndividumB = obj.matRad_bit2de(P(i,y1/2+1:y1))
                %valueIndividumB = obj.matRad_bit2de(IndividumB)
                %alpha=mod(obj.matRad_bit2de(P(i,1:y1/2)), 360)
                %alpha=mod(valueIndividumA, 360)
                %beta=mod(obj.matRad_bit2de(P(i,y1/2+1:y1)), 360)
                %beta=mod(valueIndividumB, 360)
                pln.propStf.gantryAngles
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
            x = size(P);
            Y1=zeros(p,1);
            %F = F +10; % adding 10 to ensure no chromosome has negative fitness
            % elite selection
            e=3;
            for i = 1:3
                obj.Min_fitness_value = min(F);
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

