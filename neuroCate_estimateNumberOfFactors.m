function [nZ, S, S0, ER] = neuroCate_estimateNumberOfFactors(Y, XM, estimator, nZMax)
% Produce an estimate of the number of factors using the Eigenvalue Ratio
% (ER) estimator proposed by Ahn and Horenstein (2013)  
%
% Y: data
% XM: matrix containing known covariates explaining the data
% estimator: type of estimator. For now, only ER is implemented, but future
% updates of the function may add alternative options
% nZMax: the maximum number of factors considered. It may have some
% influence on the estimation of nZ. By default, it is set to 10
if nargin == 1
  XM = [];
  estimator = 'ER';
  nZMax = 10;
elseif nargin == 2
  nZMax = 10;
  estimator = 'ER';
elseif nargin == 3
  nZMax = 10;
end

% number of known covariates
nXM = size(XM,2);

% if covariates provided, rotate the data
if  ~isempty(XM)
  [rotationMatrix, ~]= qr(XM);
  rotationMatrix = rotationMatrix';
  % rotate the data ans save it into itself for memory reason
  Y = rotationMatrix((nXM + 1):end,:) * Y;
end
[nObs, nVox] = size(Y);
m =  min(nObs, nVox);

% extract the eigenvalues of (Y'Y / (nVox * nSubj))
[~,S]  = eig(Y*Y');
S = sort(diag(S),'descend') / (nObs * nVox);

% add mock eigenvalue S0 allowing an estimation of nZ as 0; we used here
% the formula used by Ahn and Horenstein (2013) in their simulations
S0 = sum(S) / log(m);

if estimator == 'ER'
  % compute eigenvalue ratios
  ER = [S0; S(1:nZMax)] ./ S(1:(nZMax+1)) ;
end

nZ = find(ER == max(ER)) - 1;


