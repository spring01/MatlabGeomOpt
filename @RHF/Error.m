function [error, gradient] = Error(obj, coeffs, densVecSubSet, densVecRef, elecEnergyRef)
nbf = size(obj.overlapMat, 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = densVecSubSet * coeffs;

oeiVec = reshape(obj.coreHamilt, [], 1);
fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);

[densVec, elecEnergy, ~] ...
    = obj.DiagonalizeFock(reshape(fockVec, nbf, []), ...
    inv_S_Half);

% may apply weights in the future
error = norm(densVec - densVecRef) + 0*abs(elecEnergy - elecEnergyRef);

if(nargout > 1)
    gradient = zeros(length(coeffs), 1);
    delta = 1e-5;
    for i = 1:length(gradient)
        newCoeffs = coeffs;
        newCoeffs(i) = newCoeffs(i) + delta;
        gradient(i) = (obj.Error(newCoeffs, densVecSubSet, densVecRef, elecEnergyRef) - error) / delta;
    end
end
end
