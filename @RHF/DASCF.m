function [hfEnergy, iter] = DASCF(obj, iniDensVec)
nbf = size(obj.overlapMat, 1);
if(nargin < 2)
    iniDensVec = zeros(nbf^2, 1);
end
oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = iniDensVec;
elecEnergy = 0;

% diis adiis
cdiis = CDIIS(obj.overlapMat);
adiis = ADIIS(oeiVec);

cdiis_cell{1} = cdiis;
cdiis_cell{2} = CDIIS(obj.overlapMat);
cdiis_cell{3} = CDIIS(obj.overlapMat);
cdiis_cell{4} = CDIIS(obj.overlapMat);
cdiis_cell{5} = CDIIS(obj.overlapMat);
cdiis_cell{6} = CDIIS(obj.overlapMat);

for iter = 1:obj.maxSCFIter
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
        
    % diis extrapolate Fock matrix
    cdiis.Push(fockVec, densVec); % density must be idempotent
    adiis.Push(fockVec, densVec); % Fock must be built from idempotent density
    
    cdiis_cell{2}.Push(fockVec, densVec);
    cdiis_cell{3}.Push(fockVec, densVec);
    cdiis_cell{4}.Push(fockVec, densVec);
    cdiis_cell{5}.Push(fockVec, densVec);
    cdiis_cell{6}.Push(fockVec, densVec);
    if(cdiis.IAmBetter())
        densVec = cdiis.ExtrapolateDensity();
        
        % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_exp_{i} -> rho_exp_{i+1} => rho_tr_i+1
        % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_tr_{i} -> rho_exp_{i+1} => rho_tr2_i+1
        % linear regression: rho_exp_{i-4}, ... , rho_tr_{i-1}, rho_tr2_{i} -> rho_exp_{i+1} => rho_tr3_i+1
        for daggerIter = 1:length(cdiis_cell)
            cdiis_cell{daggerIter}.Push(fockVec, densVec);
            densVec = cdiis_cell{daggerIter}.ExtrapolateDensity();
        end
    else
        densVec = adiis.InterpolateDensity();
    end
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
    [densVec, elecEnergy, orbital] ...
        = obj.DiagonalizeFock(reshape(fockVec, nbf, []), ...
        inv_S_Half);
    elecEnergy = oeiVec'*densVec + elecEnergy;
    
    if(sqrt(mean((densVec - oldDensVec).^2)) < obj.RMSDensityThreshold ...
            && max(abs(densVec - oldDensVec)) < obj.MaxDensityThreshold ...
            && abs(elecEnergy - oldElecEnergy) < obj.EnergyThreshold)
        break;
    end
end
hfEnergy = elecEnergy + obj.nucRepEnergy;

obj.orbital = orbital;
obj.densVec = densVec;
obj.finalFockVec = fockVec;

end

