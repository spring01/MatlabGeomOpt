function finalCoeffs = TrainDAgger2(obj, densVecSet, refDensVecSet)

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
% maxLevel = size(refDensVecSet, 2) - 1;
maxLevel = 4;
for level = 1:maxLevel
    trainedLevels{level}.coeffs = {};
    trainedLevels{level}.densVecSet = zeros(size(densVecSet(:, 1:1+level)));
end

for iter = 2:size(refDensVecSet, 2)
    
    % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_exp_{i} -> rho_exp_{i+1} => rho_tr_i+1
    % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_tr_{i} -> rho_exp_{i+1} => rho_tr2_i+1
    % linear regression: rho_exp_{i-4}, ... , rho_tr_{i-1}, rho_tr2_{i} -> rho_exp_{i+1} => rho_tr3_i+1
    
    useIndices = iter-4:iter;
    useIndices = useIndices(useIndices > 0);
    densVecSubset = densVecSet(:,useIndices);
    refDensVec = refDensVecSet(:,iter);
    
    [trainedLevels{1}.coeffs{iter}, trainedLevels{1}.densVecSet(:, iter)] = ConstrLinReg(densVecSubset, refDensVec);
%     [tr1Coeffs{iter}, tr1DensVecSet(:, iter)] = ConstrLinReg(densVecSubset, refDensVec);
    
    newDensVecSubset = densVecSubset;
    for daggerLevel = 2:min(iter - 1, maxLevel)
        newDensVecSubset(:, end - daggerLevel + 2) = trainedLevels{daggerLevel - 1}.densVecSet(:, iter - daggerLevel + 2);
        [trainedLevels{daggerLevel}.coeffs{iter}, trainedLevels{daggerLevel}.densVecSet(:, iter)] = ConstrLinReg(newDensVecSubset, refDensVec);
    end
    
    
end

finalCoeffs = trainedLevels;

end

function [diisCoeffs, predDensVector] = ConstrLinReg(densVecSubset, refDensVec)
onesVec = ones(size(densVecSubset, 2), 1);
hessian = [densVecSubset'*densVecSubset, onesVec; onesVec', 0];
diisCoeffs = hessian \ [densVecSubset'*refDensVec; 1];
predDensVector = densVecSubset * diisCoeffs(1:end-1);
diisCoeffs = diisCoeffs(1:end-1);

% disp(norm(predDensVector - refDensVec));
% disp(norm(densVecSubset(:,end) - refDensVec));
end

