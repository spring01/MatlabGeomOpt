function [res, grad] = Target(coeffs, H)
nVecs = length(coeffs);

res = reshape(H, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;
res = reshape(res, [], nVecs) * coeffs;

if(nargout > 1)
    
    
    grad = reshape(H + permute(H, [2 1 3 4]) + permute(H, [3 2 1 4]) + permute(H, [4 2 3 1]), [], nVecs) * coeffs;
%     grad2 = reshape(4*H, [], nVecs) * coeffs;
    grad = reshape(grad, [], nVecs) * coeffs;
    grad = reshape(grad, [], nVecs) * coeffs;
end

end