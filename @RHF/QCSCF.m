function [hfEnergy, iter] = QCSCF(obj, iniDensVec)
nbf = size(obj.overlapMat, 1);
if(nargin < 2)
    iniDensVec = zeros(nbf^2, 1);
end
oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = iniDensVec;
elecEnergy = 0;

allTEIs = obj.matpsi2.Integrals_AllTEIs();

fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);

[densVec, elecEnergy, orbital] ...
    = obj.DiagonalizeFock(reshape(fockVec, nbf, []), ...
    inv_S_Half);

for iter = 1:obj.maxSCFIter
    
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
    
    
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    fockInMO = orbital' * reshape(fockVec, nbf, []) * orbital;
    
    elecEnergy = (oeiVec + fockVec)'*densVec;
    
    numOcc = obj.numElectrons/2;
    numVir = nbf - numOcc;
    
    fVec = zeros(numOcc*numVir, 1);
    for iOcc = 1:numOcc
        for aVir = 1:numVir
            fVec((iOcc-1)*numVir + aVir) = fockInMO(iOcc, numOcc+aVir);
        end
    end
    
    allTEIsInMO = TransformTensor4(allTEIs, orbital);
    AMat = zeros(numOcc*numVir);
    BMat = zeros(numOcc*numVir);
    for iOcc = 1:numOcc
        for aVir = 1:numVir
            for jOcc = 1:numOcc
                for bVir = 1:numVir
                    AMat((iOcc-1)*numVir + aVir, (jOcc-1)*numVir + bVir) = (iOcc==jOcc)*fockInMO(numOcc+aVir, numOcc+bVir) - (aVir==bVir)*fockInMO(iOcc,jOcc) + 2*allTEIsInMO(iOcc, numOcc+aVir, numOcc+bVir, jOcc) - allTEIsInMO(iOcc, jOcc, numOcc+bVir, numOcc+aVir);
                    BMat((iOcc-1)*numVir + aVir, (jOcc-1)*numVir + bVir) = 2*allTEIsInMO(iOcc, numOcc+aVir, jOcc, numOcc+bVir) - allTEIsInMO(iOcc, numOcc+bVir, jOcc, numOcc+aVir);
                end
            end
        end
    end
    [eigVec, eigVal] = eig([elecEnergy, sqrt(2)*fVec'; sqrt(2)*fVec, elecEnergy*eye(size(AMat)) + AMat + BMat]);
    diagEigVal = diag(eigVal);
    [~, minIndex] = min(diagEigVal);
    dVec = eigVec(:, minIndex);
    
    dMat = zeros(numOcc, numVir);
    for iOcc = 1:numOcc
        for aVir = 1:numVir
            dMat(iOcc, aVir) = dVec((iOcc-1)*numVir + aVir + 1);
        end
    end
    
    occOrbs = orbital(:, 1:numOcc);
    virOrbs = orbital(:, numOcc+1:end);
    for iOcc = 1:numOcc
        for aVir = 1:numVir
            occOrbs(:, iOcc) = occOrbs(:, iOcc) + dMat(iOcc, aVir) / dVec(1) * virOrbs(:, aVir);
        end
    end
    
    densVec = reshape(occOrbs * occOrbs', [], 1);
    
    orbital = occOrbs;
    for aVir = 1:numVir
        currVirOrb = virOrbs(:, aVir);
        for iOrb = 1:size(orbital, 2)
            currVirOrb = currVirOrb - (currVirOrb'*obj.overlapMat*orbital(:, iOrb)) .* orbital(:, iOrb);
        end
        orbital(:, numOcc+aVir) = currVirOrb ./ sqrt(currVirOrb'*obj.overlapMat*currVirOrb);
    end
    
%     oldDensVec = densVec;
%     oldElecEnergy = elecEnergy;
%     fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
%     [densVec, elecEnergy, orbital] ...
%         = obj.DiagonalizeFock(reshape(fockVec, nbf, []), ...
%         inv_S_Half);
%     elecEnergy = oeiVec'*densVec + elecEnergy;
    
    
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




