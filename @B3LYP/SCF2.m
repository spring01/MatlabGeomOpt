function [hfEnergy, iter] = SCF2(obj, iniDensVec)
nbf = size(obj.overlapMat, 1);
if(nargin < 2)
    iniDensVec = zeros(nbf^2, 1);
end
oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = iniDensVec;
elecEnergy = 0;

% diis adiis
comdiis = ComDIIS(obj.overlapMat, 4);
% cdiis = CDIIS(obj.overlapMat);

for iter = 1:obj.maxSCFIter
    oldDensVec = densVec;
    oldElecEnergy = elecEnergy;
    
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    
    % diis extrapolate Fock matrix
    comdiis.Push(fockVec, densVec); % density must be idempotent
%     cdiis.Push(fockVec, densVec);
%     densVec = comdiis.ExtrapolateDensity();
%     fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    fockVec = comdiis.Extrapolate();
    [densVec, elecEnergy, orbital] ...
        = obj.DiagonalizeFock(reshape(fockVec, nbf, []), ...
        inv_S_Half);
    elecEnergy = oeiVec'*densVec + elecEnergy - reshape(obj.currentV, 1, []) * densVec + obj.matpsi2.DFT_EnergyXC();
    
    obj.energySet(iter) = elecEnergy;
    
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




