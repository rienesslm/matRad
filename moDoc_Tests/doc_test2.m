function doc_test2()
        % This function illustrates a documentation test defined for MOdox.
        % Other than that it does absolutely nothing
        %
        % Examples:
        %   a=2;
        %   disp(a)
        %   % Expected output is prefixed by '%||' as in the following line:
        %   %|| 2
        %   %
        %   % The test continues because no interruption through whitespace,
        %   % as the previous line used a '%' comment character;
        %   % thus the 'a' variable is still in the namespace and can be 
        %   % accessed.
        %   b=3+a;
        %   disp(a+[3; 4; 5; 6])
        %   %|| [5; 6; 7; 8]
        %
        %   b=[5; 6; 7; 8]
        %   disp(b)
        %   %|| [5; 6; 7; 8]
        %
        %   xi = [100, 150, 200, 250, 300, 400, 500];
        %   yi = [2506.7, 2582.8, 2658.1, 2733.7, 2810.4, 2967.9, 3131.6];
        %   x = 2680.78;
        %   matRad_interp1_solution = matRad_interp1(xi, yi, x);
        %   interp1_solution = interp1(xi, yi, x);
        %   disp(matRad_interp1_solution)
        %   %|| interp1_solution
        %
        %
        % % tests end here because test indentation has ended

