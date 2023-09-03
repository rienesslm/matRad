clear;
load('Channels8_PS10_GS20_cP8_mN6_ThirdTry.mat');
load('Channels12_PS10_GS20_cP8_mN6_SecondTry.mat');
load('Channels16_PS10_GS15_cP8_mN6_SecondTry.mat')
GS5 = load('Channels16_PS10_GS20_cP8_mN6_FourthTry.mat');
load('Channels20_PS10_GS15_cP8_mN6.mat');
Channels20_5GS = load('Channels20_PS10_GS20_cP8_mN6_ThirdTry.mat');

B = [COB10GS.fitnessValueMatrix; GS5.COB10GS.fitnessValueMatrix];
C = [COB20.fitnessValueMatrix; Channels20_5GS.COB20GS.fitnessValueMatrix];

matrixList = {COB20GS.fitnessValueMatrix, COB12.fitnessValueMatrix, B, C};
min_values = {zeros(1,20), zeros(1,20), zeros(1,20), zeros(1,20)};
x = 1:length(min_values{1});

for i = 1:length(matrixList)
    min_values{i} = find_min_in_each_row(matrixList{i});
    
end

y_1 = min_values{1};
y_2 = min_values{2};
y_3 = min_values{3};
y_4 = min_values{4};
p = plot(x,y_1,'-s');
p = plot(x,y_2,'-r');
p = plot(x,y_3,'-c');
p = plot(x,y_4,'-m');


function min_values = find_min_in_each_row(matrix)
    % matrix: Input matrix
    
    % Get the number of rows in the matrix
    num_rows = size(matrix, 1);
    
    % Initialize an array to store the minimum values
    min_values = zeros(num_rows, 1);
    
    % Loop through each row and find the minimum value
    for row = 1:num_rows
        min_values(row) = min(matrix(row, :));
    end
end