load data.mat

% cdiis = CDIIS(overlap);
% for i = 1:5
%     cdiis.Push(fockVecs(:, i), densVecs(:, i))
% end
% cdiis.ExtrapolateDensity

nbf = sqrt(size(fockVecs, 1));
nVecs = size(fockVecs, 2);

S_Half = sqrtm(overlap);
e = cell(nVecs, nVecs);
for i = 1:nVecs
    for j = 1:nVecs
        fock_i = reshape(fockVecs(:, i), nbf, []);
        dens_j = reshape(densVecs(:, j), nbf, []);
        FtDt = S_Half \ fock_i * dens_j * S_Half;
        e{i, j} = reshape(FtDt - FtDt', [], 1);
    end
end

H = zeros(nVecs, nVecs, nVecs, nVecs);


for i = 1:nVecs
    for j = 1:nVecs
        for k = 1:nVecs
            for l = 1:nVecs
                H(i,j,k,l) = e{i, j}' * e{k, l};
            end
        end
    end
end


iniCoeffs = [0 0 0 0 1]';
% iniCoeffs = iniCoeffs ./ norm(iniCoeffs);

% options = optimoptions(@fmincon, 'Display', 'iter', 'GradObj', 'on');
% finalCoeffs = fmincon(@(coeffs)Target(coeffs, H), ...
%     iniCoeffs, [], [], ones(1, length(iniCoeffs)), 1, [], [], [], options);
% 
% finalCoeffs



coeffs = iniCoeffs;

% coeffs = (rand(5,1));
% coeffs = coeffs ./ sum(coeffs);
% [val, grad, hess] = Target(coeffs, H);
% onesVec = ones(length(coeffs), 1);
% hessL = [hess, onesVec; onesVec', 0];
% eig(hessL)

for iter = 1:10
[val, grad, hess] = Target(coeffs, H);
lambda = -4 * val;
gradL = [grad + lambda; 0];
if(norm(gradL) < eps)
    break;
end
onesVec = ones(length(coeffs), 1);
hessL = [hess, onesVec; onesVec', 0];
coeffsAndLambda = [coeffs; lambda];
coeffsAndLambda = coeffsAndLambda - hessL \ gradL;
coeffs = coeffsAndLambda(1:end-1);
disp(norm(gradL));
end

% 
% % temp = reshape(H + permute(H, [2 1 3 4]) + permute(H, [3 2 1 4]) + permute(H, [4 2 3 1]), [], nVecs) * coeffs;
% % temp = coeffs' * reshape(temp, nVecs, []);
% 
% for iter = 1:100
% D = zeros(nVecs, nVecs);
% for i = 1:nVecs
%     for j = 1:nVecs
%         D(i, j) = 0;
%         for k = 1:nVecs
%             for l = 1:nVecs
%                 D(i, j) = D(i, j) + coeffs(k) * coeffs(l) * (H(i,j,k,l) + H(j,i,k,l) + H(j,k,i,l) + H(j,k,l,i));
%             end
%         end
%     end
% end
% 
% lambda = coeffs' * D * coeffs;
% coeffs = D \ (lambda .* ones(nVecs, 1));
% end


