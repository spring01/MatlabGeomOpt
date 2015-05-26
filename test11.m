% cart = [16                -1.55096010   -0.68685376    0.00000000
% 1                 -1.06097692   -2.07278899    0.00000000
% 1                 -1.06095162    0.00610450    1.20025020
% 1                 -1.06095162    0.00610450   -1.20025020
% 1                 -5.55096010   -0.68680447    0.00000000];
% 
% mol = Molecule(cart);

load geomOptMols.mat
mol = mol1;
basisSet = '6-311g*';

matpsi = MatPsi2(mol.cartesian, basisSet);
scf = B3LYP.MatPsi2Interface(matpsi);
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
