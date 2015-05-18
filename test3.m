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

cell = rhfTrain.TrainDAgger2(densVecSet, refDensVecSet)


