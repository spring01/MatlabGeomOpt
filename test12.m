cart = [6                 -2.05658021    1.64108723    0.00000900
    6                 -0.65483021    1.64108723    0.00000900
    6                 -1.38788921    3.80333623   -0.00009500
    7                 -2.49630721    2.96159523    0.00000000
    1                 -3.43663521    3.25229623   -0.00021700
    1                 -2.75816821    0.80785423    0.00010700
    1                 -0.01045721    0.76838123    0.00003400
    1                 -1.49681121    4.88710823   -0.00024200
    7                 -0.23083221    3.01188323    0.00000900
    48                 1.76198223    3.70527388   -0.00003963];

mol = Molecule(cart);

basisSet = '3-21gWithCd';

matpsi = MatPsi2(mol.cartesian, basisSet);
scf = RHF.MatPsi2Interface(matpsi);
% [ener, iter] = scf.SCF();
% finalDensVec = scf.densVec;
% disp([ener, iter]);
% set1 = scf.energySet;

[ener2, iter2] = scf.SCF2();
disp([ener2, iter2]);
set2 = scf.energySet;

% plot(log10(abs(set1(end) - set1)));
% hold;
% plot(log10(abs(set2(end) - set2)), 'r');