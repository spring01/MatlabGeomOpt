function [value, grad, hess] = Target(coeffs, tensor)
nVecs = length(coeffs);

value = reshape(tensor, [], nVecs) * coeffs;
value = reshape(value, [], nVecs) * coeffs;
value = reshape(value, [], nVecs) * coeffs;
value = reshape(value, [], nVecs) * coeffs;

if(nargout > 1)
    grad = reshape( ...
        permute(tensor, [1 2 3 4]) + ...
        permute(tensor, [2 1 3 4]) + ...
        permute(tensor, [3 1 2 4]) + ...
        permute(tensor, [4 1 2 3]), [], nVecs) * coeffs;
    grad = reshape(grad, [], nVecs) * coeffs;
    grad = reshape(grad, [], nVecs) * coeffs;
end

if(nargout > 2)
    hess = reshape( ...
        permute(tensor, [1 2 3 4]) + permute(tensor, [1 3 2 4]) + permute(tensor, [1 4 2 3]) + ...
        permute(tensor, [2 1 3 4]) + permute(tensor, [2 3 1 4]) + permute(tensor, [2 4 1 3]) + ...
        permute(tensor, [3 1 2 4]) + permute(tensor, [3 2 1 4]) + permute(tensor, [3 4 1 2]) + ...
        permute(tensor, [4 1 2 3]) + permute(tensor, [4 2 1 3]) + permute(tensor, [4 3 1 2]), [], nVecs) * coeffs;
    hess = reshape(hess, [], nVecs) * coeffs;
    hess = reshape(hess, [], nVecs);
end

end