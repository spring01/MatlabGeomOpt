load geomOptMols.mat

basisSet = '6-31g*';

molIni = mol1;
matpsiIni = MatPsi2(molIni.cartesian, basisSet);
matpsiIni.SCF_RunSCF();
iniDensVec = reshape(matpsiIni.SCF_DensityAlpha(), [], 1);

molTrain = mol2;
matpsiTrain = MatPsi2(molTrain.cartesian, basisSet);
rhfTrain = RHF.MatPsi2Interface(matpsiTrain);
[ener, iter] = rhfTrain.SCF(iniDensVec);
finalDensVec = rhfTrain.densVec;
disp([ener, iter]);

[ener2, iter2] = rhfTrain.SCF2(iniDensVec);
disp([ener2, iter2]);
