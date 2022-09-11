function tests = test_DoseObjectives
test_functions = {};
MyFolderInfo = dir('..\optimization\+DoseObjectives');
for i=1:length(MyFolderInfo)
  if (MyFolderInfo(i).name(1)~='.')
    if not(isequal(MyFolderInfo(i).name, 'matRad_DoseObjective.m'))
      test_functions{end+1} = @() testDoseGrad(MyFolderInfo(i).name);
      test_functions{end+1} = @() testDoseObjectiveFunction(MyFolderInfo(i).name);
    end
  end
end
test_functions = test_functions';
tests = functiontests(test_functions);
end
    
function testDoseGrad(testCase, fileName)
    concatenationOfFileName = strcat('DoseObjectives.', fileName);
    functionNameHandle = str2func(concatenationOfFileName(1:end-2));
    dose = [10 20 30 40 50 60]';
    obj = functionNameHandle(); %Constructor
    doseGrad=obj.computeDoseObjectiveGradient(dose);
    verifyTrue(testCase, isreal(doseGrad)); %nth root problem (Problem based on complex values of EUD)
    [grad,err,finaldelta]=gradest(@(x) obj.computeDoseObjectiveFunction(x'), dose');
    epsilon = doseGrad*obj.maxDerivativeError + err';
    %assertElementsAlmostEqual(grad', doseGrad, 'absolute', max(epsilon)); % Accepts vectors and matrices | Accepts vectors as epsilon
    %assertVectorsAlmostEqual(grad', doseGrad, 'absolute', max(epsilon)); % Just accepts verctors | Add err + errImp
    verifyEqual(testCase, size(doseGrad), size(dose));
    verifyNotEqual(testCase, doseGrad, NaN);
    verifyNotEqual(testCase, doseGrad, Inf);
    verifyNotEqual(testCase, doseGrad, -Inf);
end

function testDoseObjectiveFunction(testCase, fileName)
    concatenationOfFileName = strcat('DoseObjectives.', fileName);
    functionNameHandle = str2func(concatenationOfFileName(1:end-2));
    dose = [10 20 30 40 50 60];
    obj = functionNameHandle(); %Constructor
    f = obj.computeDoseObjectiveFunction(dose);
    [r,c] = size(f);
    verifyTrue(testCase, isreal(f));
    verifyGreaterThan(testCase, f, zeros(r,c) - eps);
    verifyNotEqual(testCase, f, NaN);
    verifyNotEqual(testCase, f, Inf);
    verifyNotEqual(testCase, f, -Inf);
end