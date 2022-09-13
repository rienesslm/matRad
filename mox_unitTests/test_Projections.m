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
    test_functions = {};
    MyFolderInfo = dir('../optimization/projections'); %relative path
    for i=1:length(MyFolderInfo)
      if (MyFolderInfo(i).name(1)~='.')
       if not(isequal(MyFolderInfo(i).name, 'matRad_BackProjection.m'))
        test_functions{end+1} = eval( ...
            ['@() testProjection(''' MyFolderInfo(i).name ''', dij, w )']);
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

    %Fails for matRad_VariableRBEProjection..
    g = proj.projectSingleScenarioGradient(dij,{ones(nVox,1)},1,w);
    
    [jacobEst,err] = jacobianest(@(x) proj.computeSingleScenario(dij,1,x'),w');
    gEst = sum(jacobEst)'; %Sums each row of the matrix(jacobEst)
    err = sum(err); % 1x2 vector
    assertEqual(size(g), size(w));
    assertEqual(size(g), size(gEst));
    assertEqual(length(g), nB);
    assertElementsAlmostEqual(g, gEst, 'absolute', max(err));
    assertTrue(isreal(g));
    assertNotEqual(g, NaN);
    assertNotEqual(g, Inf);
    assertNotEqual(g, -Inf);
