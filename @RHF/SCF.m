function [hfEnergy, iter] = SCF(obj, iniDensMat)
nbf = size(obj.overlapMat, 1);
if(nargin < 2)
    iniDensMat = zeros(nbf, nbf);
end
oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = reshape(iniDensMat, [], 1);
elecEnergy = 0;

% diis adiis
cdiis = CDIIS(obj.overlapMat);
adiis = ADIIS(oeiVec);

iniGMat = 2 .* obj.jkFactory.JK_DensToJ(iniDensMat) ... % +2J
    - obj.jkFactory.JK_DensToK(iniDensMat); % -K

fockVec = oeiVec + reshape(iniGMat, [], 1);
fockSimVec = fockVec;
% obj.jkFactory.JK_Initialize('PKJK');
for iter = 1:obj.maxSCFIter
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
    [densVec, elecEnergy] ...
        = DiagonalizeFock(reshape(fockSimVec, nbf, []), ...
        inv_S_Half, obj.numElectrons);
    elecEnergy = oeiVec'*densVec + elecEnergy;
    
    if(sqrt(mean((densVec - oldDensVec).^2)) < obj.RMSDensityThreshold ...
            && max(abs(densVec - oldDensVec)) < obj.MaxDensityThreshold ...
            && abs(elecEnergy - oldElecEnergy) < obj.EnergyThreshold)
        break;
    end
    densMat = reshape(densVec, nbf, []);
    gMat = 2 .* obj.jkFactory.JK_DensToJ(densMat) ... % +2J
        - obj.jkFactory.JK_DensToK(densMat); % -K
    fockVec = oeiVec + reshape(gMat, [], 1);
        
    % diis extropolate Fock matrix
    cdiis.Push(fockVec, densVec); % density must be idempotent
    adiis.Push(fockVec, densVec); % Fock must be built from idempotent density
    if(cdiis.IAmBetter())
        fockSimVec = cdiis.Extrapolate();
    else
        fockSimVec = adiis.Interpolate();
    end
end
hfEnergy = elecEnergy + obj.nucRepEnergy;

obj.hfEnergy = hfEnergy;
obj.densMat = reshape(densVec, nbf, []);

end


function [densityVec, elecEnergy] = DiagonalizeFock(fockMat, inv_S_Half, numElectrons)
[orbitalOtho, orbitalEnergies] = eig(inv_S_Half*fockMat*inv_S_Half);
[orbitalEnergies, ascend_order] = sort(diag(orbitalEnergies));
orbital = inv_S_Half * orbitalOtho(:, ascend_order);
densityVec = reshape( ...
    orbital(:, 1:numElectrons/2) * orbital(:, 1:numElectrons/2)', ...
    [], 1);
elecEnergy = sum(orbitalEnergies(1:numElectrons/2));
end

