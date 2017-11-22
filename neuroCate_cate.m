function [results] = neuroCate_cate(Y,X,M,nZ,varargin)
%   'inference' - Define the type of inference
%       Choices are:
%           'parametric' - only a parametric inference is conducted (the default).
%           'permutation'- A permutation inference is conducted (the parametric
%                          inference comes for free)
%           'wb'         - A Wild Bootstrap inference is conducted (the parametric
%                          inference comes for free)
%   'saveResampledFScores>=' - Per default, does not save anything
%        (i.e. behave like the set value is Inf)
%        In practice, those scores can be used to compute a non-parameric
%        FDR as described in Millstein and Volfson (2013). For example,
%        setting this to 4, would save only F-scores greater or equal to 4.
%        The importance of not saving everything is to limit the size of the
%        output.
%   'saveResampledPValues<=' - Per default, does not save anything
%        (i.e. behave like the set value is 0)
%        In practice, those p-values can be used to compute a non-parameric
%        FDR as described in Millstein and Volfson (2013). For example,
%        setting this to 0.01 would save all resampled uncorrected p-values 
%        smaller than or equal to 0.01. 
%        The importance of not saving everything is to limit the size of the
%        output. 
%   'clusterFormingPValueThreshold' - Per default, does not do any
%        cluster-wise inference. If a value is set and a non-parametric
%        procedure is used, a non-parametric cluster-wise analysis will be
%        run.
%        Note that a p-value threshold is expexted here
%    'XYZ_vox' - Location of the voxels needed to detect clusters for 
%        cluster-wise inference. 

% deal with the optional parameters
paramNames = {'inference', 'residAdj', 'nB', 'seed', 'resamplingMatrix', 'saveResampledFScores>=', 'saveResampledPValues<=', 'clusterFormingPValueThreshold', 'XYZ_vox'};
defaults   = {'parametric', true, 999, 0, [], [], [], [], []};

[inference, residAdj, nB, seed, resamplingMatrix, minScoreFToSave, maxPValueToSave, pClusThresh, XYZ_vox]...
  = internal.stats.parseArgs(paramNames, defaults, varargin{:});

% compute useful variables from inputs
nX = size(X,2);
nM = size(M,2);
if (nM > 1)
  error('Further checks needed to validate neuroCate with several covariates of interests')
end
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
    if strcmpi(inference, 'permutation') || strcmpi(inference, 'perm')
      resamplingMatrix = zeros(nSubj, nB);
      for iB = 1:nB
        resamplingMatrix(:,iB) = randperm(nSubj);
      end
    elseif strcmpi(inference, 'wild bootstrap') || strcmpi(inference, 'wb')
      resamplingMatrix = randi([0 1], nSubj, nB)*2 - 1;
    end
    
  else % need to check if the supplied resamling matrix looks ok
    if ~all(size(resamplingMatrix) == [nSubj, nB])
      error('The specified resampling matrix does not have the appropriate dimensions. Please revise your analysis specification.')
    end
    if strcmpi(inference, 'permutation') || strcmpi(inference, 'perm')
      for iB = 1:nB
        if any(sort(resamplingMatrix(:,iB)) ~= (1:nSubj)')
          error('The specified resampling matrix is not an appropriate permutation resampling matrix');
        end
      end
    elseif strcmpi(inference, 'wild bootstrap') || strcmpi(inference, 'wb')
      for iB = 1:nB
        if any(sort(unique(resamplingMatrix(:,iB))) ~= [-1;1])
          error('The specified resampling matrix is not an appropriate Wild Bootstrap resampling matrix');
        end
      end
    end
    if isempty(pClusThresh) && isempty(XYZ_vox)
        clusterWise = false;
    elseif isempty(pClusThresh) || isempty(XYZ_vox)
        error('There is not enough information to conduct a cluster-wise inference. Please set both the voxel locations and cluster-forming threshold')
    else
        clusterWise = true;
        if pClusThresh > 0.5
            error('The cluster-wise forming threshold does not seem small enough. As a reminder, this should be a p-value threshold. We recommend a default value of 0.001')
        end
        % convert p-value threshold to F treshold for efficiency
        fClusThresh = tinv(pClusThresh/2, approxDof)^2;
    end
  end
else
    clusterWise = false;
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
  [rotZZ, betaZ, sigmaSquare] = neuroCate_fastPCA(rotY((nM + 1):end,:), nZ);
  
  % run the robust regression
  scaledRotMY = rotatedMM \ rotY(1:nM,:); % scale the data and variance estimate by (Q_M' M)^-1
  scaleForSigmaSquare = diag(inv(rotatedMM' * rotatedMM));%TODO in R cate package, they do diag(inv(rotatedMM * rotatedMM'))????
  alphaM = neuroCate_robustRegression(scaledRotMY, betaZ,...
    sigmaSquare, scaleForSigmaSquare, 100, 10^-10);
  clear scaledRotY
  
  % estimate betaM
  betaM = scaledRotMY - alphaM * betaZ;
  if nM == 1
    % varBetaM = ((OmegaM + alphaM * alphaM') / nSubj) * sigmaSquare ;
    scoreF = (nSubj/ (OmegaM + alphaM * alphaM')) * (betaM.^2 ./ sigmaSquare);
  else
    % TODO
  end
else % naive regression
  Z = [];
  betaZ = [];
  sigmaSquare = sum((rotY((nM + 1):end,:)).^2, 1) / approxDof;
  alphaM =[];
  betaM = rotatedMM \ rotY(1:nM,:);
  if nM == 1
    % varBetaM = ((OmegaM + alphaM * alphaM') / nSubj) * sigmaSquare ;
    scoreF = (nSubj/OmegaM) * (betaM.^2 ./ sigmaSquare);
  else
    % TODO
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute parametric p-values
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nM == 1
  %pParamUnc = fcdf(scoreF, 1, approxDof, 'upper');
  pParamUnc = tcdf(-sqrt(scoreF), approxDof)*2; % useful for older version of matlab
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
if strcmpi(inference, 'permutation') || strcmpi(inference, 'perm') || strcmpi(inference, 'wild bootsrap') || strcmpi(inference, 'wb')
  maxScoreF = zeros(nB+1,1);
  maxScoreF(1) = max(scoreF);
  pNonParamUnc = ones(1, nVox);
  pNonParamFWE = ones(1, nVox);
  
  % save cluster sizes for cluster-wise inferences
  if clusterWise
    maxClusterSize = nan(nB+1,1);
    activatedVoxels = scoreF >= fClusThresh;
    LocActivatedVoxels = XYZ_vox(:,activatedVoxels);
    if isempty(LocActivatedVoxels)
      maxClusterSize(1) = 0;
      nCluster          = 0;
      clusterSize       = [];
      pClusNonParamFWE  = [];
      clusterAssignment = [];
      warning('No voxels survived the cluster forming thresholding.');
    else
      clusterAssignment = spm_clusters(LocActivatedVoxels);
      nCluster     = max(clusterAssignment);
      clusterSize = histc(clusterAssignment,1:nCluster);
      maxClusterSize(1) = max(clusterSize);
      pClusNonParamFWE = ones(1, nCluster);
    end
  end
    
  % create cell arrays to save resampled score or p-values if required
  if ~isempty(minScoreFToSave)
  	results.resampledScores = cell(nB, 1);
  end
  if ~isempty(maxPValueToSave)
    results.resampledPValues = cell(nB, 1);
  end
  
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
    % some code to print the progress
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
      if strcmpi(inference, 'permutation') || strcmpi(inference, 'perm')
        rotY_b = fittedRotY + ...
          (rotationMatrix((nX+1):end,:) * rotationMatrix((nX+1):end, resamplingMatrix(:,iB))') * rotResid;
      elseif strcmpi(inference, 'wild bootstrap') || strcmpi(inference, 'wb')
        rotY_b = fittedRotY + ...
          ((rotationMatrix((nX+1):end,:) .* kron(resamplingMatrix(:,iB)', ones(nSubj-nX,1))) * rotationMatrix((nX+1):end,:)') * rotResid;
      else
        error('unknown inference procedure');
      end
    else % naive regression
      if strcmpi(inference, 'permutation') || strcmpi(inference, 'perm')
        rotY_b = ...
          (rotationMatrix((nX+1):end,:) * rotationMatrix((nX+1):end, resamplingMatrix(:,iB))') * rotY;
      elseif strcmpi(inference, 'wild bootstrap') || strcmpi(inference, 'wb')
        rotY_b = ...
          ((rotationMatrix((nX+1):end,:) .* kron(resamplingMatrix(:,iB)', ones(nSubj-nX,1))) * rotationMatrix((nX+1):end,:)') * rotY;
      else
        error('unknown inference procedure');
      end
    end
    
    % estimate the model an compute resampled scores
    if nZ > 0
      % run the factor analysis using PCA
      [~, betaZ_b, sigmaSquare_b] = neuroCate_fastPCA(rotY_b((nM + 1):end,:), nZ);
      
      % run the robust regression
      scaledRotMY_b = rotatedMM \ rotY_b(1:nM,:); % scale the data and variance estimate by (Q_M' M)^-1
      alphaM_b = neuroCate_robustRegression(scaledRotMY_b, betaZ_b,...
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
        
    % Save some resampled scores or p-values if required
    if ~isempty(minScoreFToSave)
        results.resampledScores{iB} = scoreF_b(scoreF_b >= minScoreFToSave);
    end
    if ~isempty(maxPValueToSave)
        % first, need to convert to p-values 
        p_b = tcdf(-sqrt(scoreF_b), approxDof)*2;
        results.resampledPValues{iB} = p_b(p_b <= maxPValueToSave);
    end
    if clusterWise
      activatedVoxels_b = scoreF_b >= fClusThresh;
      LocActivatedVoxels_b = XYZ_vox(:,activatedVoxels_b);
      if isempty(LocActivatedVoxels_b)
        maxClusterSize(iB + 1) = 0;
      else
        clusterAssignment_b = spm_clusters(LocActivatedVoxels_b);
        nCluster_b     = max(clusterAssignment_b);
        clusterSize_b = histc(clusterAssignment_b,1:nCluster_b);
        maxClusterSize(iB + 1) = max(clusterSize_b);
      end
      if ~isempty(LocActivatedVoxels)
        pClusNonParamFWE = pClusNonParamFWE + (maxClusterSize(iB + 1) >= clusterSize);
      end
    end   
  end % for iB = 1:nB
  
  % format p-values correctly
  pNonParamUnc = pNonParamUnc / (nB + 1);
  pNonParamFWE = pNonParamFWE / (nB + 1);
  if clusterWise
    results.pClusNonParamFWE  = pClusNonParamFWE / (nB + 1);
    results.maxClusterSize    = maxClusterSize;
    results.clusterSize       = clusterSize;
    results.clusterAssignment = nan(1, nVox);
    results.clusterAssignment(activatedVoxels) = clusterAssignment;
  end
  % save results
  results.pNonParamUnc = pNonParamUnc;
  results.pNonParamFWE = pNonParamFWE;
  results.maxScoreF    = maxScoreF;
  fprintf('\n');
end


