function tests = test_matRad_interp1
tests = functiontests(localfunctions);
end

function testInterp1Extrapolation(testCase)
    xi = [100, 150, 200, 250, 300, 400, 500];
    yi = [2506.7, 2582.8, 2658.1, 2733.7, 2810.4, 2967.9, 3131.6]';
    x = 650;
    matRad_interp1_solution = matRad_interp1(xi, yi, x);
    interp1_solution = interp1(xi, yi, x);
    verifyEqual(testCase, matRad_interp1_solution, interp1_solution)
end

function testInterp1Interpolation(testCase)
    xi = [100, 150, 200, 250, 300, 400, 500];
    yi = [2506.7, 2582.8, 2658.1, 2733.7, 2810.4, 2967.9, 3131.6]';
    x = 222;
    matRad_interp1_solution = matRad_interp1(xi, yi, x);
    interp1_solution = interp1(xi, yi, x);
    verifyEqual(testCase, matRad_interp1_solution,interp1_solution)
end
