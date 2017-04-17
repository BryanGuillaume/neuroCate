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
  
   % check some inputs regarding non-paramteric procedures before going further
  if ~strcmpi(inference, 'parametric') && ~strcmpi(inference, 'permutation') && ~strcmpi(inference, 'wb')
    error('Unknown inference procedure. Please revise your analysis specification.')
  end
  if ~strcmpi(inference, 'parametric') 
    if isempty(resamplingMatrix) % build it from the seed
    	% set the random seed for repeatablity
      rng(seed);
      % build the resampling matrix now
      if strcmpi(inference, 'permutation')
        resamplingMatrix = zeros(nSubj, nB);
        for iB = 1:nB
          resamplingMatrix(:,iB) = randperm(nSubj);
        end
      else % Wild Bootstrap resampling matrix from the Rademacher distr.
        resamplingMatrix = randi([0 1], nSubj, nB)*2 - 1;
      end
      
    else % need to check if the supplied resamling matrix looks ok
      if ~all(size(resamplingMatrix) == [nSubj, nB])
        error('The specified resampling matrix does not have the appropriate dimensions. Please revise your analysis specification.')       
      end
      if strcmpi(inference, 'permutation')
        for iB = 1:nB
          if any(sort(resamplingMatrix(:,iB)) ~= (1:nSubj)')
            error('The specified resampling matrix is not an appropriate permutation resampling matrix');
          end
        end      
      else % WB resampling matrix
        for iB = 1:nB
          if any(sort(unique(resamplingMatrix(:,iB))) ~= [-1;1])
            error('The specified resampling matrix is not an appropriate Wild Bootstrap resampling matrix');
          end
        end     
      end
    end
  end
  
  % compute the rotation matrix and related variables
  [rotationMatrix, ~]= qr(XM);
  rotationMatrix = rotationMatrix';
  rotatedMM = rotationMatrix((nX + 1):(nX + nM),:) * M; % Q_M' M
  rotY = rotationMatrix((nX + 1):end,:) * Y;
  
  % remove Y for memory
  clear Y
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % estimation of the orginal model parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  if nZ > 0
    % run the factor analysis using PCA
    [rotZZ, betaZ, sigmaSquare] = fastPCA(rotY((nM + 1):end,:), nZ);
    
    % run the robust regression
    scaledRotMY = rotatedMM \ rotY(1:nM,:); % scale the data and variance estimate by (Q_M' M)^-1    
    scaleForSigmaSquare = diag(inv(rotatedMM' * rotatedMM));%TODO in R cate package, they do diag(inv(rotatedMM * rotatedMM'))????
    alphaM = robustRegression(scaledRotMY, betaZ,...
      sigmaSquare, scaleForSigmaSquare, 100, 10^-10);
    clear scaledRotY
    
    % estimate betaM
    betaM = scaledRotMY - alphaM * betaZ;

   else % naive regression
    Z = [];
    betaZ = [];
    sigmaSquare = sum((rotY((nM + 1):end,:)).^2, 1) / approxDof;
    alphaM =[];
    betaM = rotatedMM \ rotY(1:nM,:);
  end
  
  
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
        rotResid = rotResid ./ kron(sqrt(1-diag(rotMZZ * pinv(rotMZZ))), ones(1,nVox));
      end
    else % if nZ = 0, deal with this later
%       fittedRotY <- 0 * rotY
%       rotResid <- rotY
    end
    fprintf('Bootsrap # 0'); 
    for iB = 1:nB
      if iB <= 10
        fprintf('\b%i', iB); 
      elseif iB <= 100
        fprintf('\b\b%i', iB); 
      elseif iB <= 1000
        fprintf('\b\b\b%i', iB); 
      elseif iB <= 10000
        fprintf('\b\b\b\b%i', iB);
      end

      % resample the data
      if nZ > 0
        if strcmpi(inference, 'permutation')
          rotY_b = fittedRotY + ...
            (rotationMatrix((nX+1):end,:) * rotationMatrix((nX+1):end, resamplingMatrix(:,iB))') * rotResid;
        elseif strcmpi(inference, 'wb')
          rotY_b = fittedRotY + ...
            ((rotationMatrix((nX+1):end,:) .* kron(resamplingMatrix(:,iB)', ones(nSubj-nX,1))) * rotationMatrix((nX+1):end,:)') * rotResid;
        else
          error('unknown inference procedure');
        end
      else % naive regression
        if strcmpi(inference, 'permutation')
          rotY_b = ...
            (rotationMatrix * rotationMatrix((nX+1):end, resamplingMatrix(:,iB))') * rotY;
        elseif strcmpi(inference, 'wb')
          rotY_b = ...
            ((rotationMatrix .* kron(resamplingMatrix(:,iB), ones(1,nSubj))) * rotationMatrix((nX+1):end,:)') * rotY;
        else
          error('unknown inference procedure');
        end
      end
      
      % estimate the model an compute resampled scores
      if nZ > 0
        % run the factor analysis using PCA
        [~, betaZ_b, sigmaSquare_b] = fastPCA(rotY_b((nM + 1):end,:), nZ);

        % run the robust regression
        scaledRotMY_b = rotatedMM \ rotY_b(1:nM,:); % scale the data and variance estimate by (Q_M' M)^-1    
        alphaM_b = robustRegression(scaledRotMY_b, betaZ_b,...
          sigmaSquare_b, scaleForSigmaSquare, 100, 10^-10);

        % estimate betaM
        betaM_b = scaledRotMY_b - alphaM_b * betaZ_b;
        if nM == 1
          scoreF_b = (nSubj/ (OmegaM + alphaM_b * alphaM_b')) * (betaM_b.^2 ./ sigmaSquare_b);
        else
          % TODO when nM > 1, do some vectorising trick here    
        end
      else % naive regression
        sigmaSquare_b = sum((rotY_b((nM + 1):end,:)).^2, 1) / approxDof;
        betaM_b = rotatedMM \ rotY_b(1:nM,:);
        if nM == 1
          scoreF_b = (nSubj/OmegaM) * (betaM_b.^2 ./ sigmaSquare_b);
        else
          % TODO when nM > 1, do some vectorising trick here    
        end
      end

      % update maxScore and non-parametric p-values
      maxScoreF(iB + 1) = max(scoreF_b);
      pNonParamUnc = pNonParamUnc + (scoreF_b >= scoreF);
      pNonParamFWE = pNonParamFWE + (maxScoreF(iB + 1) >= scoreF);
    end % for iB = 1:nB
    
    % format p-values correctly
    pNonParamUnc = pNonParamUnc / (nB + 1);
    pNonParamFWE = pNonParamFWE / (nB + 1);
    
    % save results
    results.pNonParamUnc = pNonParamUnc;
    results.pNonParamFWE = pNonParamFWE;
    results.maxScoreF    = maxScoreF;
    fprintf('\n'); 
  end
end

