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
        e{i, j} = FtDt - FtDt';
    end
end

H = zeros(nVecs, nVecs, nVecs, nVecs);


for i = 1:nVecs
    for j = 1:nVecs
        for k = 1:nVecs
            for l = 1:nVecs
                H(i,j,k,l) = H(i,j,k,l) + sum(sum(e{i, j} .* e{k, l}));
            end
        end
    end
end


iniCoeffs = [0 0 0 0 1]';
% iniCoeffs = iniCoeffs ./ norm(iniCoeffs);

options = optimoptions(@fmincon, 'Display', 'iter', 'GradObj', 'off');
finalCoeffs = fmincon(@(coeffs)Target(coeffs, H), ...
    iniCoeffs, [], [], ones(1, length(iniCoeffs)), 1, [], [], [], options);

finalCoeffs


