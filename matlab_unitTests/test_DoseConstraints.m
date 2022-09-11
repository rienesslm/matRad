function tests = test_DoseConstraints
test_functions = {};
MyFolderInfo = dir('..\optimization\+DoseConstraints'); %relative path
for i=1:length(MyFolderInfo)
  if (MyFolderInfo(i).name(1)~='.')
    if not(isequal(MyFolderInfo(i).name, 'matRad_DoseConstraint.m'))
      test_functions{end+1} = eval(['@() testComputeDoseConstraintJacobian(''' MyFolderInfo(i).name ''')']); %Too show which object failed. Too have object name expliceed in function name.
      test_functions{end+1} = eval(['@() testComputeDoseConstraintFunction(''' MyFolderInfo(i).name ''')']);
    end
  end
end
test_functions = test_functions';
tests = functiontests(test_functions);
end
    
function testComputeDoseConstraintJacobian(testCase, fileName)
    concatenationOfFileName = strcat('DoseConstraints.', fileName);
    functionNameHandle = str2func(concatenationOfFileName(1:end-2));
    dose = [10 20 30 40 50 60]';
    obj = functionNameHandle(); %Constructor
    doseConstraintJacobian=obj.computeDoseConstraintJacobian(dose);
    [jac,err] = jacobianest(@(x) obj.computeDoseConstraintFunction(x'), dose');
    epsilon = obj.maxDerivativeError + max(err);
    verifyTrue(testCase, isreal(doseConstraintJacobian)); %nth root
    %assertVectorsAlmostEqual(jac', doseConstraintJacobian, 'relative', epsilon); % Build 'assertVectorsAlmostEqual' through subtraction 'LessThan'
    verifyEqual(testCase, size(doseConstraintJacobian), size(dose));
    verifyNotEqual(testCase, doseConstraintJacobian, NaN);
    verifyNotEqual(testCase, doseConstraintJacobian, Inf);
    verifyNotEqual(testCase, doseConstraintJacobian, -Inf);
end

function testComputeDoseConstraintFunction(fileName)
    concatenationOfFileName = strcat('DoseConstraints.', fileName);
    functionNameHandle = str2func(concatenationOfFileName(1:end-2));
    dose = [10 20 30 40 50 60];
    obj = functionNameHandle(); %Constructor
    f = obj.computeDoseConstraintFunction(dose);
    [r,c] = size(f);
    verifyTrue(testCase, isreal(f));
    verifyGreaterThan(testCase, f, zeros(r,c) - eps);
    verifyNotEqual(testCase, f, NaN);
    verifyNotEqual(testCase, f, Inf);
    verifyNotEqual(testCase, f, -Inf);
end