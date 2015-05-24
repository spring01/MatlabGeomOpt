cart = [16                -1.55096010   -0.68685376    0.00000000
1                 -1.06097692   -2.07278899    0.00000000
1                 -1.06095162    0.00610450    1.20025020
1                 -1.06095162    0.00610450   -1.20025020
1                 -5.55096010   -0.68680447    0.00000000];

mol = Molecule(cart);

basisSet = 'cc-pvdz';

molTrain = mol;
matpsiTrain = MatPsi2(molTrain.cartesian, basisSet);
rhfTrain = RHF.MatPsi2Interface(matpsiTrain);
[ener, iter] = rhfTrain.SCF();
finalDensVec = rhfTrain.densVec;
disp([ener, iter]);

[ener2, iter2] = rhfTrain.SCF2();
disp([ener2, iter2]);

