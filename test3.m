load geomOptMolsBig.mat

basisSet = 'sto-3g';

molIni = mol1;
matpsiIni = MatPsi2(molIni.cartesian, basisSet);
matpsiIni.SCF_RunSCF();
iniDensVec = reshape(matpsiIni.SCF_DensityAlpha(), [], 1);


molTrain = mol2;
matpsiTrain = MatPsi2(molTrain.cartesian, basisSet);

rhfTrain = RHF.MatPsi2Interface(matpsiTrain);
[ener, densVecSet, refDensVecSet, iter] = rhfTrain.ExpertSCF(iniDensVec);
disp([ener, iter]);

trainedLevels = rhfTrain.TrainDAgger2(densVecSet, refDensVecSet);

[ener2, iter2] = rhfTrain.TestDAggerSCF(trainedLevels, iniDensVec);
disp([ener2, iter2]);



molTest = mol3;
matpsiTest = MatPsi2(molTest.cartesian, basisSet);
rhfTest = RHF.MatPsi2Interface(matpsiTest);
[ener3, iter3] = rhfTest.SCF(iniDensVec);
disp([ener3, iter3]);
[ener4, iter4] = rhfTest.TestDAggerSCF(trainedLevels, iniDensVec);
disp([ener4, iter4]);

