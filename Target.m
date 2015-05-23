function [res, grad, grad2] = Target(coeffs, H)
nVecs = length(coeffs);

res = reshape(H, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;

if(nargout > 1)
    grad = zeros(size(coeffs));
    for i = 1:length(coeffs)
        delta = zeros(size(coeffs));
        delta(i) = 1e-7;
        grad(i) = (Target(coeffs + delta, H) - res) ./ delta(i);
    end
    
    
    grad2 = reshape(H + permute(H, [2 1 3 4]) + permute(H, [3 2 1 4]) + permute(H, [4 2 3 1]), [], nVecs) * coeffs;
%     grad2 = reshape(4*H, [], nVecs) * coeffs;
    grad2 = reshape(grad2, [], nVecs) * coeffs;
    grad2 = reshape(grad2, [], nVecs) * coeffs;
end

end