function [hfEnergy, coeffs, iter] = TrainDAggerSCF(obj, iniDensVec, collectionDensVec)

% keep only the densities we essentially use
stepSize = 2;
finalDensVec = collectionDensVec(:,end);
numDensVectors = size(collectionDensVec, 2);
collectionDensVec = collectionDensVec(:,1:stepSize:end);
if(mod(numDensVectors, 2))
    collectionDensVec = [collectionDensVec, finalDensVec];
end

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

for iter = 1:size(collectionDensVec, 2)
    
    % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_exp_{i} -> rho_exp_{i+1} => rho_tr_i+1
    % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_tr_{i} -> rho_exp_{i+1} => rho_tr2_i+1
    % linear regression: rho_exp_{i-4}, ... , rho_tr_{i-1}, rho_tr2_{i} -> rho_exp_{i+1} => rho_tr3_i+1
    for daggerIter = 1:iter-1
    end
    
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    
    
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
    [densVec, elecEnergy, orbital] ...
        = DiagonalizeFock(reshape(fockVec, nbf, []), ...
        inv_S_Half, obj.numElectrons);
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

end


function [densityVec, elecEnergy, orbital, orbEigValues] ...
    = DiagonalizeFock(fockMat, inv_S_Half, numElectrons)
[orbitalOtho, orbEigValues] = eig(inv_S_Half*fockMat*inv_S_Half);
[orbEigValues, ascend_order] = sort(diag(orbEigValues));
orbital = inv_S_Half * orbitalOtho(:, ascend_order);
densityVec = reshape( ...
    orbital(:, 1:numElectrons/2) * orbital(:, 1:numElectrons/2)', ...
    [], 1);
elecEnergy = sum(orbEigValues(1:numElectrons/2));
end

