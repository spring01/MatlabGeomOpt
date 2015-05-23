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


molTest = mol3;
matpsiTest = MatPsi2(molTest.cartesian, basisSet);
rhfTest = RHF.MatPsi2Interface(matpsiTest);

randMat = rand(sqrt(numel(finalDensVec))) - 0.5;
randMat = randMat + randMat';
randVec = 0.1 .* reshape(randMat, [], 1);
[ener3, iter3] = rhfTrain.SCF(finalDensVec + randVec);
disp([ener3, iter3]);
