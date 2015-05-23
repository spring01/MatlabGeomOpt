load sih4.mat

basisSet = '6-31+g*';

molTrain = mol;
matpsiTrain = MatPsi2(molTrain.cartesian, basisSet);
rhfTrain = RHF.MatPsi2Interface(matpsiTrain);
[ener, iter] = rhfTrain.SCF();
finalDensVec = rhfTrain.densVec;
disp([ener, iter]);

[ener2, iter2] = rhfTrain.SCF2();
disp([ener2, iter2]);