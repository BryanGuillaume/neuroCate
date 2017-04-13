function [results] = cate(Y,X,M,nZ,varargin)
  %    'inference'  - Define the type of inference
  %    'parametric' - only a parametric inference is conducted (the default).
  %    'permutation'- A permutation inference is conducted (the parametric 
  %                   inference comes for free)
  %    'wb'         - A Wild Bootstrap inference is conducted (the parametric 
  %                   inference comes for free) 
  
  % deal with the optional parameters
  paramNames = {'inference', 'residAdj', 'nB', 'seed', 'resamplingMatrix'};
  defaults   = {'parametric', true, 999, 0, []};

  [inference, residAdj, nB, seed, resamplingMatrix]...
      = internal.stats.parseArgs(paramNames, defaults, varargin{:});  
    
  % compute useful variables from inputs
  nX = size(X,2);
  nM = size(M,2);
  [nSubj, nVox] = size(Y);
  XM = [X M];
  invSigmaXM = inv(XM' * XM / nSubj);
  OmegaM = invSigmaXM((nX + 1):end, (nX + 1):end);
  approxDof = nSubj - nX - nM - nZ;
  
  % compute the rotation matrix and related variables
  [rotationMatrix, ~]= qr(XM);
  rotationMatrix = rotationMatrix';
  rotatedMM = rotationMatrix((nX + 1):(nX + nM),:) * M; % Q_M' M
  rotY = rotationMatrix((nX + 1):end,:) * Y;
  
  % remove Y for memory
  clear Y
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % run the factor analysis using PCA
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  if nZ > 0
%     [U,S,V] = svd(rotY((nM + 1):end,:),'econ');
%     rotZZ =  sqrt(nSubj - nX - nM) * U(:,1:nZ);
%     betaZ = S(1:nZ,1:nZ) * V(:,1:nZ)' / sqrt(nSubj - nX - nM);
%     sigmaSquare = sum((rotY((nM + 1):end,:) - rotZZ * betaZ).^2, 1) / approxDof;
%     clear U S V
%     
    [rotZZ, betaZ, sigmaSquare] = fastPCA(rotY((nM + 1):end,:), nZ);
    
  else % naive regression
    Z = [];
    betaZ = [];
    sigmaSquare = sum((rotY((nM + 1):end,:)).^2, 1) / approxDof;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % run the robust regression to estimate alphaM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % scale the data and variance estimate by (Q_M' M)^-1
  scaledRotMY = rotatedMM \ rotY(1:nM,:);
  %TODO in R cate package, they do diag(inv(rotatedMM * rotatedMM'))????
  scaleForSigmaSquare = diag(inv(rotatedMM' * rotatedMM));
  
  % run a robust regression for each m = 1,...,nM
%   alphaM = zeros(nM,nZ);
%   kHuber = 1.345;
%   sqrtWeight = zeros(1, nVox);
%   for iM = 1:nM
%     scaledSigma = sqrt(sigmaSquare * scaleForSigmaSquare(iM));
%     convergence = false;
%     convergenceThreshold = 10^-10;
%     nIter = 100;
%     % starting point using an OLS estimator
%     alphaM(iM,:) = scaledRotMY(iM,:) * pinv(betaZ')';
%     iter = 1;
%     while(~convergence && iter <= nIter)
%       alphaM_previous = alphaM(iM,:);
%       scaledResid = (scaledRotMY(iM,:) - alphaM(iM,:) * betaZ) ./ scaledSigma;
%       selection = (abs(scaledResid) <= kHuber);
%       sqrtWeight(selection) = 1;
%       sqrtWeight(~selection) = sqrt(kHuber ./ abs(scaledResid(~selection)));
%       
%       % weight variable using the sqrt of the Huber weight
%       wScaledRotMY = scaledRotMY(iM,:) .* sqrtWeight;
%       wBetaZ = betaZ .* sqrtWeight(ones(nZ,1),:);      
%       alphaM(iM,:) = wScaledRotMY(iM,:) * pinv(wBetaZ')';
% 
%       % check convergence
%       convergence = max(abs(1 - alphaM(iM,:)./alphaM_previous)) < convergenceThreshold;
%       iter = iter + 1;
%     end
%     if ~convergence
%       warning('The robust regression has failed to converge.')
%     end     
%   end % end of the robust regression
  alphaM = robustRegression(scaledRotMY, betaZ,...
    sigmaSquare, scaleForSigmaSquare, 100, 10^-10);
  clear wBetaZ wScaledRotMY sqrtWeight scaledSigma

  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % estimate betaM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  betaM = scaledRotMY - alphaM * betaZ;
  
  % clear large variables for memory
  clear scaledRotY
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % compute scores
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % assumes identity contrast matrix 
  
  if nM == 1
    % varBetaM = ((OmegaM + alphaM * alphaM') / nSubj) * sigmaSquare ;
    scoreF = (nSubj/ (OmegaM + alphaM * alphaM')) * (betaM.^2 ./ sigmaSquare);
  else
    
    % TODO when nM > 1, do some vectorising trick here
    
    % scoreF = (OmegaM + alphaM * alphaM')
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % compute parametric p-values
  %%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  if nM == 1
    pParamUnc = fcdf(scoreF, 1, approxDof, 'upper');
  else
    % TODO when nM > 1
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % compute Z
  %%%%%%%%%%%%%%%%%%%%%%%%%%%   
  if nZ > 0
    Z = rotationMatrix' * [rotationMatrix(1:(nX+nM),:) * M * alphaM; rotZZ];
  end
  
  results.pParamUnc = pParamUnc;
  results.scoreF = scoreF;
  results.betaM = betaM;
  results.alphaM = alphaM;  
  results.betaZ = betaZ;
  results.sigmaSquare = sigmaSquare;
  results.Z = Z;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % non-parametric inference
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmpi(inference, 'permutation') || strcmpi(inference, 'wb')
    maxScoreF = zeros(nB+1,1);
    maxScoreF(1) = max(scoreF);
    pNonParamUnc = ones(1, nVox);
    pNonParamFWE = ones(1, nVox);
    
    % compute the fitted data + residuals of the original data
    if nZ > 0
      rotMZZ = [rotatedMM * alphaM ; rotZZ];
      fittedRotY = rotMZZ * betaZ;
      rotResid = rotY - fittedRotY;
      clear rotY
      if residAdj
        rotResid = rotResid ./ kron(sqrt(1-diag(rotMZZ * pinv(rotMZZ))), ones(1,10));
      end
    else % if nZ = 0, deal with this later
%       fittedRotY <- 0 * rotY
%       rotResid <- rotY
    end  
    for iB = 1:nB
      
      if strcmpi(inference, 'permutation')
    
      elseif strcmpi(inference, 'wb')
    
      end
      
      betaM_b = scaledRotMY_b - alphaM_b * betaZ_b;
      
      if nM == 1
      % varBetaM = ((OmegaM + alphaM * alphaM') / nSubj) * sigmaSquare ;
        scoreF_b = (nSubj/ (OmegaM + alphaM_b * alphaM_b')) * (betaM_b.^2 ./ sigmaSquare_b);
      else
        % TODO when nM > 1, do some vectorising trick here    
      end
      maxScoreF(iB + 1) = max(scoreF_b)
    end
      
  end
 
end

