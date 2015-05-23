function res = Target(coeffs, H)
nVecs = length(coeffs);

res = reshape(H, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;

end