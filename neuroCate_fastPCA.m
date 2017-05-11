function [rotZZ, betaZ, sigmaSquare] = neuroCate_fastPCA(Y, nZ)
% function for solving the factor models in cate
% Instead being based on an SVD on Y or an eigenvalue decompostion of Y'*Y,
% the fast PCA is based on a truncated eigenvalue decompostion of Y*Y'
nRow = size(Y,1);
[UZ,~] = eigs(Y * Y', nZ);
rotZZ = sqrt(nRow) * UZ;
betaZ = (UZ' * Y) / sqrt(nRow);
sigmaSquare = sum((Y- rotZZ * betaZ).^2, 1) / (nRow - nZ);

