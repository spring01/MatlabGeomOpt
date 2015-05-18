function finalCoeffs = TrainDAgger(obj, densVecSet, elecEnergySet)

% keep only the densities we essentially use
stepSize = 2;
if(mod(size(densVecSet, 2), 2) == 0)
    densVecSet = [densVecSet, densVecSet(:,end)];
    elecEnergySet = [elecEnergySet, elecEnergySet(end)];
end
densVecSet = densVecSet(:,1:stepSize:end);
elecEnergySet = elecEnergySet(:,1:stepSize:end);

finalCoeffs = {};
for iter = 1:size(densVecSet, 2)-1
    
    % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_exp_{i} -> rho_exp_{i+1} => rho_tr_i+1
    % linear regression: rho_exp_{i-4}, ... , rho_exp_{i-1}, rho_tr_{i} -> rho_exp_{i+1} => rho_tr2_i+1
    % linear regression: rho_exp_{i-4}, ... , rho_tr_{i-1}, rho_tr2_{i} -> rho_exp_{i+1} => rho_tr3_i+1
    
    useIndices = iter-8:iter;
    useIndices = useIndices(useIndices > 0);
    densVecSubSet = densVecSet(:,useIndices);
    densVecRef = densVecSet(:,iter+1);
    elecEnergyRef = elecEnergySet(iter+1);
    
    iniCoeffs = zeros(length(useIndices), 1);
    iniCoeffs(end) = 1;
    
    options = optimoptions(@fmincon, 'Display', 'iter', 'GradObj', 'off');
    finalCoeffs{iter} = fmincon(@(coeffs)obj.Error(coeffs, densVecSubSet, densVecRef, elecEnergyRef), ...
        iniCoeffs, ones(1, length(iniCoeffs)), 0, ones(1, length(iniCoeffs)), 1, [], [], [], options);
    
    
    for daggerIter = 1:iter-1
    end
    
    
end

end

