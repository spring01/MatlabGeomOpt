function trainedLevels = TrainDAgger2(obj, densVecSet, refDensVecSet)
nbf = size(obj.overlapMat, 1);
oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

% keep only the densities we essentially use
stepSize = 2;
densVecSet = [densVecSet(:, 1), densVecSet];
% if(mod(size(densVecSet, 2), 2) == 0)
%     densVecSet = [densVecSet, densVecSet(:,end)];
% end
densVecSet = densVecSet(:, 1:stepSize:end);

% if(mod(size(refDensVecSet, 2), 2) == 0)
%     refDensVecSet = [refDensVecSet, refDensVecSet(:,end)];
% end
refDensVecSet = refDensVecSet(:, 1:stepSize:end);


trainedLevels = {};
maxLevel = size(refDensVecSet, 2) - 2;
% maxLevel = 7;
for level = 1:maxLevel
    trainedLevels{level}.coeffs = {};
    trainedLevels{level}.densVecSet = zeros(size(densVecSet, 1), 0);
end

for iter = 2:size(refDensVecSet, 2)
    
    % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_exp_{i} -> rho_exp_{i+1} => rho_tr_i+1
    % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_tr_{i} -> rho_exp_{i+1} => rho_tr2_i+1
    % linear regression: rho_exp_{i-4}, ... , rho_tr_{i-1}, rho_tr2_{i} -> rho_exp_{i+1} => rho_tr3_i+1
    
    useIndices = iter-8:iter;
    useIndices = useIndices(useIndices > 0);
    densVecSubset = densVecSet(:,useIndices);
    refDensVec = refDensVecSet(:,iter);
    
    
    for daggerLevel = 1:min(iter - 1, maxLevel)
        newDensVecSubset = densVecSubset;
        for replaceLevel = 1:daggerLevel-1
            newDensVecSubset(:, end + replaceLevel - daggerLevel + 1) = trainedLevels{replaceLevel}.densVecSet(:, iter + replaceLevel - daggerLevel + 1);
        end
        [trainedLevels{daggerLevel}.coeffs{iter}, densVecSim] = ConstrLinReg(newDensVecSubset, refDensVec);
        fockVecSim = oeiVec + reshape(obj.DensToG(reshape(densVecSim, nbf, [])), [], 1);
        [solvedDensVec, ~, ~] ...
            = obj.DiagonalizeFock(reshape(fockVecSim, nbf, []), ...
            inv_S_Half);
        trainedLevels{daggerLevel}.densVecSet(:, iter+1) = solvedDensVec;
    end
    
%     disp('one iter')
end

end

function [diisCoeffs, predDensVector] = ConstrLinReg(densVecSubset, refDensVec)
onesVec = ones(size(densVecSubset, 2), 1);
hessian = [densVecSubset'*densVecSubset, onesVec; onesVec', 0];
diisCoeffs = hessian \ [densVecSubset'*refDensVec; 1];
predDensVector = densVecSubset * diisCoeffs(1:end-1);
diisCoeffs = diisCoeffs(1:end-1);

% disp(norm(predDensVector - refDensVec));
% disp(norm(densVecSubset(:,end) - refDensVec));
% disp('one level')
end

