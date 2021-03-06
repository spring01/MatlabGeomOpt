function [hfEnergy, densVecSet, refDensVecSet, iter] = ExpertSCF(obj, iniDensVec)
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

densVecSet = densVec;
refDensVecSet = zeros(length(densVec), 0);

for iter = 1:obj.maxSCFIter
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    
    % diis extrapolate Fock matrix
    cdiis.Push(fockVec, densVec); % density must be idempotent
    densVec = cdiis.ExtrapolateDensity();
    refDensVecSet(:, iter) = densVec;
    
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    
    [densVec, elecEnergy, orbital] ...
        = obj.DiagonalizeFock(reshape(fockVec, nbf, []), ...
        inv_S_Half);
    elecEnergy = oeiVec'*densVec + elecEnergy;
    
    % collection of densities
    densVecSet(:, iter+1) = densVec;
    
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

