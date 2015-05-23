function tensor4 = TransformTensor4(tensor4, trans)
nbf1 = size(trans, 1);
nbf2 = size(trans,2);

% pqrs -> iqrs
tensor4 = reshape(trans' * reshape(tensor4, nbf1,[]), nbf2,nbf1,nbf1,nbf1);

% iqrs -> ijrs
tensor4 = permute(tensor4,[2 1 3 4]);
tensor4 = trans' * reshape(tensor4, nbf1,[]);
tensor4 = permute(reshape(tensor4, nbf2,nbf2,nbf1,nbf1), [2 1 3 4]);

% ijrs -> ijks
tensor4 = permute(tensor4,[1 2 4 3]);
tensor4 = reshape(tensor4, [], nbf1) * trans;
tensor4 = permute(reshape(tensor4, nbf2,nbf2,nbf1,nbf2), [1 2 4 3]);

% ijks -> ijkl
tensor4 = reshape(tensor4, [], nbf1) * trans;
tensor4 = reshape(tensor4, nbf2,nbf2,nbf2,nbf2);
end