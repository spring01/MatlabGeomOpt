load mols.mat

mol = mols{22};
matpsi = MatPsi2(mol.cartesian, '6-31g*');

molTest = mols{23};
matpsiTest = MatPsi2(molTest.cartesian, '6-31g*');

molIni = mols{21};
matpsiIni = MatPsi2(molIni.cartesian, '6-31g*');
matpsiIni.SCF_RunSCF();
iniDensVec = reshape(matpsiIni.SCF_DensityAlpha(), [], 1);

rhf = RHF.MatPsi2Interface(matpsi);
[ener, iter] = rhf.SCF(iniDensVec);
disp([ener, iter]);

finalDensVec = rhf.densVec;

% [enerExpert, collection, iterExpert] = rhf.ExpertSCF();

[ener2, allCoeffs, iter2] = rhf.CheatSCF(finalDensVec, iniDensVec);
disp([ener2, iter2]);


rhfTest = RHF.MatPsi2Interface(matpsiTest);

[ener3, iter3] = rhfTest.TestCheatSCF(allCoeffs, iniDensVec);
disp([ener3, iter3]);
