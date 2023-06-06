classdef matRad_GeneticOptimizer
    %MATRAD_GENETICOPTIMIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties     
        population=100; % population size
        chromosomePairs=30; % number of pairs of chromosomes to be crossovered
        mutationNumber=30; % number chromosomes to be mutated
        numberGenerations=250; % Total number of generations
    end
    
    methods
        function obj = matRad_GeneticOptimizer()
            figure
            title('Blue - Average         Red-Maximum');
            xlabel('Genration')
            ylabel('Objective Function Value')
            hold on
            P=obj.initPopulation(obj.population);
            K=0;
            [x1,y1]=size(P);
            P1=0;
            for i=1:obj.numberGenerations
                Cr=obj.crossover(P, obj.chromosomePairs);
                Mu=obj.mutation(P, obj.mutationNumber);
                P(obj.population + 1 : obj.population + 2 * obj.chromosomePairs, :)=Cr;
                P(obj.population + 2 * obj.chromosomePairs + 1 : obj.population + 2 * obj.chromosomePairs + obj.mutationNumber, :)=Mu;
                E=evaluation(obj, P);
                [P,S]=selection(obj, P, E, obj.population); %No clue why the function call needs 'obj' here...
                K(i,1)=sum(S)/obj.population;
                K(i,2)=S(1); %best
                plot(K(:,1),'b.'); drawnow
                hold on
                plot(K(:,2),'r.'); drawnow
            end
            Max_fitness_value=max(K(:,2))
            P2 = P(1,:); % Best chromosome
            % convert binary to real number
            
            A=obj.matRad_bit2de(P2(1, 1:y1/2));
            x=-3+A*(3-(-3))/(2^(y1/2)-1);
            B=obj.matRad_bit2de(P2(1, y1/2+1:y1));
            y=-3+B*(3-(-3))/(2^(y1/2)-1);
            Optimal_solution=[x,y]
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
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
            Y=round(rand(n,40));
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
                A1(1,r2:y1)=A2(1,r2:40);
                A2(1,r2:40)=B1;
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

        function Y=evaluation(obj, P)
            [x1, y1]=size(P);
            H=zeros(1,x1);
            for i = 1:x1
                A=obj.matRad_bit2de(P(i,1:y1/2));
                x=-3+A*(3-(-3))/(2^(y1/2)-1);
                B=obj.matRad_bit2de(P(i,y1/2+1:y1));
                y=-3+B*(3-(-3))/(2^(y1/2)-1);
                % Maximization problem
                H(1,i)=3*(1-x)^2*exp(-x^2 - (y+1)^2)...
                    - 10*(x/5 - x^3 - y^5)*exp(-x^2 - y^2)...
                    -1/3*exp(-(x+1)^2 - y^2);
            end
            Y=H;
        end

        function [YY1, YY2] = selection(obj, P, F, p)
            % P - population, F -fitness value, p - population size
            [x,y]=size(P);
            Y1=zeros(p,y);
            F = F +10; % adding 10 to ensure no chromosome has negative fitness
            % elite selection
            e=3;
            for i = 1:3
                [r1,c1]=find(F==max(F));
                Y1(i,:)=P(max(c1),:);
                P(max(c1),:)=[];
                Fn(i)=F(max(c1));
                F(:,max(c1))=[];
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
            YY2=Fn-10; % substract 10 to return the original fitness
        end
        
        function decimalNumber = matRad_bit2de(obj, binaryString)
            % Convert binary strings to decimal numbers
            decimalNumber = polyval(binaryString, 2);
        end
    end
end

