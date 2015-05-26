load geomOptMols.mat

basisSet = '6-31g*';

molIni = mol1;
matpsiIni = MatPsi2(molIni.cartesian, basisSet);
matpsiIni.SCF_RunSCF();
iniDensVec = reshape(matpsiIni.SCF_DensityAlpha(), [], 1);

molTrain = mol2;
matpsi = MatPsi2(molTrain.cartesian, basisSet);
scf = RHF.MatPsi2Interface(matpsi);
[ener, iter] = scf.SCF();
finalDensVec = scf.densVec;
disp([ener, iter]);
set1 = scf.energySet;

[ener2, iter2] = scf.SCF2();
disp([ener2, iter2]);
set2 = scf.energySet;

plot(log10(abs(set1(end) - set1)));
hold;
plot(log10(abs(set2(end) - set2)), 'r');

