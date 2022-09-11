function test_suite=test_DoseConstraints
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
    initTestSuite;
    

function testComputeDoseConstraintJacobian(fileName)
    concatenationOfFileName = strcat('DoseConstraints.', fileName);
    functionNameHandle = str2func(concatenationOfFileName(1:end-2));
    dose = [10 20 30 40 50 60]';
    obj = functionNameHandle(); %Constructor
    doseConstraintJacobian=obj.computeDoseConstraintJacobian(dose);
    [jac,err] = jacobianest(@(x) obj.computeDoseConstraintFunction(x'), dose');
    epsilon = obj.maxDerivativeError + max(err);
    assertTrue(isreal(doseConstraintJacobian)); %nth root
    assertVectorsAlmostEqual(jac', doseConstraintJacobian, 'relative', epsilon); % Just accepts verctors | Add err + errImp
    assertEqual(size(doseConstraintJacobian), size(dose));
    assertNotEqual(doseConstraintJacobian, NaN);
    assertNotEqual(doseConstraintJacobian, Inf);
    assertNotEqual(doseConstraintJacobian, -Inf);

function testComputeDoseConstraintFunction(fileName)
    concatenationOfFileName = strcat('DoseConstraints.', fileName);
    functionNameHandle = str2func(concatenationOfFileName(1:end-2));
    dose = [10 20 30 40 50 60];
    obj = functionNameHandle(); %Constructor
    f = obj.computeDoseConstraintFunction(dose);
    [r,c] = size(f);
    assertTrue(isreal(f));
    assertGreaterThan(f, zeros(r,c) - eps);
    assertNotEqual(f, NaN);
    assertNotEqual(f, Inf);
    assertNotEqual(f, -Inf);

