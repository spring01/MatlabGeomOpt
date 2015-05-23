% load mols.mat
% mol1 = mols{2};
% mol3 = mols{3};

load geomOptMols.mat

basisSet = 'sto-3g';

molIni = mol1;
matpsiIni = MatPsi2(molIni.cartesian, basisSet);
matpsiIni.SCF_RunSCF();
iniDensVec = reshape(matpsiIni.SCF_DensityAlpha(), [], 1);
% iniDensVec = iniDensVec + 0.2*rand(size(iniDensVec));

molTrain = mol3;
matpsiTrain = MatPsi2(molTrain.cartesian, basisSet);
rhfTrain = RHF.MatPsi2Interface(matpsiTrain);
[ener, iter] = rhfTrain.SCF(iniDensVec);
finalDensVec = rhfTrain.densVec;
disp([ener, iter]);

[ener2, iter2] = rhfTrain.QCSCF(iniDensVec);
disp([ener2, iter2]);
