function [hfEnergy, iter] = TestCheatSCF(obj, allCoeffs, iniDensVec)
nbf = size(obj.overlapMat, 1);
if(nargin < 3)
    iniDensVec = zeros(nbf^2, 1);
end
oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = iniDensVec;
elecEnergy = 0;

allDensVectors = zeros(nbf^2, 0);

for iter = 1:obj.maxSCFIter
    
    try
        allDensVectors(:,iter) = densVec;
        temp = zeros(nbf^2, iter+5);
        temp(:, end - size(allDensVectors, 2) + 1 : end) = allDensVectors;
        densVec = temp(:, end - length(allCoeffs{iter}) + 1 : end) * allCoeffs{iter};
    catch
        
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

