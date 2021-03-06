function [hfEnergy, iter] = TestDAggerSCF(obj, trainedLevels, iniDensVec)
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

stepCoeffs{1} = 1;
for step = 2:length(trainedLevels)
%     stepCoeffs{step} = trainedLevels{step-1}.coeffs{step};
    stepCoeffs{step} = trainedLevels{1}.coeffs{step};
end

for iter = 1:obj.maxSCFIter
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    
    % diis extrapolate Fock matrix
    cdiis.Push(fockVec, densVec); % density must be idempotent
    if(iter < length(stepCoeffs))
        densVec = densVecSet(:, end-length(stepCoeffs{iter})+1:end) * stepCoeffs{iter};
    else
        densVec = cdiis.ExtrapolateDensity();
    end
    
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
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

