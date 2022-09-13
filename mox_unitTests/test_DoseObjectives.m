function test_suite=test_DoseObjectives
    test_functions = {};
    MyFolderInfo = dir('../optimization/+DoseObjectives'); %relative path
    for i=1:length(MyFolderInfo)
      if (MyFolderInfo(i).name(1)~='.')
       if not(isequal(MyFolderInfo(i).name, 'matRad_DoseObjective.m'))
        test_functions{end+1} = eval(['@() testDoseGrad(''' MyFolderInfo(i).name ''')']); %Too show which object failed. Too have object name expliceed in function name.
        test_functions{end+1} = eval(['@() testDoseObjectiveFunction(''' MyFolderInfo(i).name ''')']);
       end
      end
    end
    test_functions = test_functions';
    initTestSuite;
    

function testDoseGrad(fileName)
    concatenationOfFileName = strcat('DoseObjectives.', fileName);
    functionNameHandle = str2func(concatenationOfFileName(1:end-2));
    dose = [10 20 30 40 50 60]';
    obj = functionNameHandle(); %Constructor
    doseGrad=obj.computeDoseObjectiveGradient(dose);
    assertTrue(isreal(doseGrad)); %nth root problem (Problem based on complex values of EUD)
    [grad,err,finaldelta]=gradest(@(x) obj.computeDoseObjectiveFunction(x'), dose');
    epsilon = doseGrad*obj.maxDerivativeError + err';
    assertElementsAlmostEqual(grad', doseGrad, 'absolute', max(epsilon)); % Accepts vectors and matrices | Accepts vectors as epsilon
    assertVectorsAlmostEqual(grad', doseGrad, 'absolute', max(epsilon)); % Just accepts verctors | Add err + errImp
    assertEqual(size(doseGrad), size(dose));
    assertNotEqual(doseGrad, NaN);
    assertNotEqual(doseGrad, Inf);
    assertNotEqual(doseGrad, -Inf);

function testDoseObjectiveFunction(fileName)
    concatenationOfFileName = strcat('DoseObjectives.', fileName);
    functionNameHandle = str2func(concatenationOfFileName(1:end-2));
    dose = [10 20 30 40 50 60];
    obj = functionNameHandle(); %Constructor
    f = obj.computeDoseObjectiveFunction(dose);
    [r,c] = size(f);
    assertTrue(isreal(f));
    assertGreaterThan(f, zeros(r,c) - eps);
    assertNotEqual(f, NaN);
    assertNotEqual(f, Inf);
    assertNotEqual(f, -Inf);