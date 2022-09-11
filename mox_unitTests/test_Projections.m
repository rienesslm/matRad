function test_suite=test_Projections
    nVox = 3;
    nB = 2;
    dij.ax= 0.1*ones(nVox,1);
    dij.bx=0.05*ones(nVox,1);
    dij.gamma=dij.ax./(2*dij.bx);
    dij.physicalDose{1} = rand(nVox,nB);
    dij.mAlphaDose{1} = rand(nVox,nB);
    dij.mLetDose{1} = rand(nVox, nB);
    dij.mSqrtBetaDose{1} = rand(nVox, nB);
    dij.ixDose = dij.bx~=0;
    dij.doseGrid.numOfVoxels = nVox;
    dij.fixedCurrent = 300;
    dij.RBE = 1.1;
    w = ones(nB,1);
    %cst{1,3}= 'TARGET';
    %cst{1,4}{1} = 1:nVox;
    test_functions = {};
    MyFolderInfo = dir('..\optimization\projections'); %relative path
    for i=1:length(MyFolderInfo)
      if (MyFolderInfo(i).name(1)~='.')
       if not(isequal(MyFolderInfo(i).name, 'matRad_BackProjection.m'))
               %isequal(MyFolderInfo(i).name, 'matRad_VariableRBEProjection.m'))
               %isequal(MyFolderInfo(i).name, 'matRad_EffectProjection.m') | ...
               %isequal(MyFolderInfo(i).name, 'matRad_ConstantRBEProjection.m') | ...
               %4isequal(MyFolderInfo(i).name, 'matRad_DoseProjection.m') | ...
               %isequal(MyFolderInfo(i).name, 'matRad_VariableRBEProjection.m'))
        test_functions{end+1} = eval(['@() testProjection(''' MyFolderInfo(i).name ''', dij, w )']); %Too show which object failed. Too have object name expliceed in function name.
        %test_functions{end+1} = @() testProjection(MyFolderInfo(i).name, dij, w);
       end
      end
    end
    test_functions = test_functions';
    initTestSuite;

function testProjection(fileName,dij,w)
    functionNameHandle = str2func(fileName(1:end-2));
    proj = functionNameHandle(); %Constructor
    nVox = dij.doseGrid.numOfVoxels;
    nB = numel(w);
    g = proj.projectSingleScenarioGradient(dij,{ones(nVox,1)},1,w); %Fails for matRad_VariableRBEProjection..
    
    name=class(proj);
    [jacobEst,err] = jacobianest(@(x) proj.computeSingleScenario(dij,1,x'),w');
    gEst = sum(jacobEst)'; %Sums each row of the matrix(jacobEst)
    err = sum(err); % 1x2 vector
    gDiff = g - gEst;
    assertEqual(size(g), size(w));
    assertEqual(size(g), size(gEst));
    assertEqual(length(g), nB);
    %assertEqual(size(jacobEst), [(length(g)+1), length(g)]); % size(jacobEst) == m x n (n==length(g), m==length(g) +1 ?)
    assertElementsAlmostEqual(g, gEst, 'absolute', max(err)); % abs(gEst.*err) == 2x1 vector * 1x2 vector...result 2x2 vector => Not comparable...
    assertTrue(isreal(g));
    assertNotEqual(g, NaN);
    assertNotEqual(g, Inf);
    assertNotEqual(g, -Inf);
