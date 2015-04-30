load mols.mat

mol = mols{2};

opt = RHFGeomOpt(mol, '6-31g*');
opt.DoGeomOpt

