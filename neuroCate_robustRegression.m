function [alphaM] = neuroCate_robustRegression(scaledRotMY, betaZ, sigmaSquare, scaleForSigmaSquare, nIter, convergenceThreshold)
% function that extimate alphaM through robust regression (Huber loss function assumed)
[nM, nVox] = size(scaledRotMY);
nZ = size(betaZ, 1);
alphaM = zeros(nM,nZ);
kHuber = 1.345;
sqrtWeight = zeros(1, nVox);
for iM = 1:nM
  scaledSigma = sqrt(sigmaSquare * scaleForSigmaSquare(iM));
  convergence = false;
  % starting point using an OLS estimator
  alphaM(iM,:) = scaledRotMY(iM,:) * pinv(betaZ')';
  iter = 1;
  while(~convergence && iter <= nIter)
    alphaM_previous = alphaM(iM,:);
    scaledResid = (scaledRotMY(iM,:) - alphaM(iM,:) * betaZ) ./ scaledSigma;
    selection = (abs(scaledResid) <= kHuber);
    sqrtWeight(selection) = 1;
    sqrtWeight(~selection) = sqrt(kHuber ./ abs(scaledResid(~selection)));
    
    % weight variable using the sqrt of the Huber weight
    wScaledRotMY = scaledRotMY(iM,:) .* sqrtWeight;
    wBetaZ = betaZ .* sqrtWeight(ones(nZ,1),:);
    alphaM(iM,:) = wScaledRotMY(iM,:) * pinv(wBetaZ')';
    
    % check convergence
    convergence = max(abs(1 - alphaM(iM,:)./alphaM_previous)) < convergenceThreshold;
    iter = iter + 1;
  end
  if ~convergence
    warning('The robust regression has failed to converge.')
  end
end % end of the robust regression