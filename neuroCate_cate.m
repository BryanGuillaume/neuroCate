function [results] = neuroCate_cate(Y,X,M,nZ,varargin)
% MANDATORY INPUTS:
%   Y - output data. It must be an nSubj by nVox matrix. This matrix can be 
%         specified as a variable in the workspace or as a ".mat" file
%         containing (only) such a matrix
%   X - Design  matrix of nX covariates of no interest. It must be an
%         nSubj by nX matrix specified as a variable in the workspace or as
%         a ".mat" file containing (only) such a matrix
%   M - Design  matrix of 1 covariate of direct interest. It must be an
%         nSubj by 1 matrix specified as a variable in the workspace or as
%         a ".mat" file containing (only) such a matrix
%   nZ - Number of unknown covariates to account for. It can be estimated
%         using the function neuroCate_estimateNumberOfFactors.m. 
% OPTIONAL INPUTS (must be specified in pair):        
%   'inference' - Define the type of inference
%         Choices are:
%           'parametric' - only a parametric inference is conducted (the default).
%           'permutation'- A permutation inference is conducted (the parametric
%                          is also computed)
%           'wb'         - A Wild Bootstrap inference is conducted (the parametric
%                          is also computed)
%   'residAdj' - boolean value (true or false) that indicates if the
%         resampled residuals are adjusted for small sample bias. Per
%         default, this is set to true.
%   'nB' - number of resamples (does not count the original sample). Per
%         defaults, it is set to 999. However, we recommend to use a higher
%         number if possible, such as 9,999 in order to increase the precision of
%         the estimated p-values. 
%   'seed' - seed number used for the resampling if a resampling matrix is
%         not supplied. Per default, it is set to 0.
%   'resamplingMatrix' - resampled matrix supplied to ensure that is
%         consistent accross several run of CATE. This is important for an EWAS,
%         particularly if setting the seed cannot be trusted (e.g., when the systems 
%         where CATE is run vary). Its size must be nSubj by nB. It can be computed
%         using the function neuroCate_resamplingMatrix.m and can be
%         specified either as a variable in the workspace or as a ".mat"
%         file containing (only) the resampling matrix.
%   'saveResampledFScores>=' - Save all the F-scores >= to the set values.
%         Per default, does not save anything (i.e. behave like the set value is Inf)
%         In practice, those F-scores can be used to compute a non-parameric
%         FDR as described in Millstein and Volfson (2013). For example,
%         setting this to 4, would save only F-scores greater or equal to 4.
%         The importance of not saving everything is to limit the size of the
%         output.
%   'saveResampledPValues<=' - Save all P-values <= to the set P-value. 
%         Per default, does not save anything (i.e. behave like the set
%         value is 0). In practice, those p-values can be used to compute 
%         a non-parameric FDR as described in Millstein and Volfson (2013). 
%         For example, setting this to 0.01 would save all resampled uncorrected  
%         p-values smaller than or equal to 0.01. 
%         The importance of not saving everything is to limit the size of the
%          output. 
%   'clusterFormingPValueThreshold' - Per default, the function does not do any
%         cluster-wise inference. If a value is set and a non-parametric
%         procedure is used, a non-parametric cluster-wise analysis will be
%         run.
%         Note that a p-value threshold is expected here.
%    'XYZ_vox' - Location of the voxels needed to detect clusters for 
%         cluster-wise inference.
%    'outputName' - Per default, the function does not save anything on disk 
%         and only return the results on the workspace. If a name is specified, 
%         the output will also be saved on disk using this name.

% deal with the optional parameters
paramNames = {'inference', 'residAdj', 'nB', 'seed', 'resamplingMatrix', 'saveResampledFScores>=', 'saveResampledPValues<=', 'clusterFormingPValueThreshold', 'XYZ_vox', 'outputName'};
defaults   = {'parametric', true, 999, 0, [], [], [], [], [], []};
try
    [inference, residAdj, nB, seed, resamplingMatrix, minScoreFToSave, maxPValueToSave, pClusThresh, XYZ_vox, outputName]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});
catch
  for i = 1:length(paramNames)
    if any(strcmpi(paramNames{i}, varargin))
      % change the default value to the one specified in varargin
      defaults{i} = varargin{(find(strcmpi(paramNames{i}, varargin)) + 1)};
    end
  end
  % parse the values in the parameters
  inference         = defaults{1};
  residAdj          = defaults{2};
  nB                = defaults{3};
  seed              = defaults{4};
  resamplingMatrix  = defaults{5};
  minScoreFToSave   = defaults{6};
  maxPValueToSave   = defaults{7};
  pClusThresh       = defaults{8};
  XYZ_vox           = defaults{9};
  outputName        = defaults{10};
end

% load inputs if specified as path to ".mat" files
if ischar(X)
  X = importdata(X);
end
if ischar(M)
  M = importdata(M);
end
if ischar(Y)
  Y = importdata(Y);
end
if ischar(resamplingMatrix)
  resamplingMatrix = importdata(resamplingMatrix);
end
if ischar(XYZ_vox)
  XYZ_vox = importdata(XYZ_vox);
end

% if some inputs are characters instead of numeric/logical  values, convert them
if ischar(residAdj)
  if strcmpi(residAdj, 'true')
    residAdj = true;
  elseif strcmpi(residAdj, 'false')
    residAdj = false;
  else
    error('residAdj is not well specified. Please use true or false')
  end
end
if ischar(nZ)
  nZ = str2double(nZ);
end
if ischar(nB)
  nB = str2double(nB);
end
if ischar(seed)
  seed = str2double(seed);
end
if ischar(minScoreFToSave)
  minScoreFToSave = str2double(minScoreFToSave);
end
if ischar(maxPValueToSave)
  maxPValueToSave = str2double(maxPValueToSave);
end
if ischar(pClusThresh)
  pClusThresh = str2double(pClusThresh);
end

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

% check if the inputs have a correct sizes
if size(X,1) ~= size(M,1) || size(X,1) ~= nSubj
  error('Inputs not well specified! Please ensure that the subjects are specified in rows.')
end

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
try
  [rotationMatrix, ~]= qr(XM);
catch
  [rotationMatrix, varToDelete]= qr(XM);
  clear varToDelete;
end
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
      try
        [~, betaZ_b, sigmaSquare_b] = neuroCate_fastPCA(rotY_b((nM + 1):end,:), nZ);
      catch
        [varToDelete, betaZ_b, sigmaSquare_b] = neuroCate_fastPCA(rotY_b((nM + 1):end,:), nZ);
        clear varToDelete;
      end
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
    results.pClusThresh = pClusThresh;
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
  
  if ~isempty(outputName)
    save(outputName, 'results');
  end
end


