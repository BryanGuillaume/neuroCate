% This script explains how neuroCate can be used to conduct an EWAS 
% analysis using CATE. It is written assuming an EWAS, but it can be 
% easily adapted for other type of mass-univariate analysis using CATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brief description of the variables:
% nSubj:  number of subjets
% nVox:   number of voxels (note that we use the term voxels here because we 
%           assume that the data comes from 3D images, but any other type of 
%           data should work; for example, vertex)
% nZ:     number of unknown covariates
% Y:      neuroimaging data matrix of size nSubj by nVox
% X:      design matrix of size nSubj by nX containing nX covariates of no
%           direct interest
% M:      methylation data matrix of size nSubj by nML containing data for 
%           nML DNA methylation loci
% Z:      design matrix of size nSubj by nZ containing nZ unknown covariates
% nB:     number of resampling (the original data is not counted in this number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Load or format the input data Y in a correct nSubj X nVox matrix format if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The piece of code below uses some SPM functions to transform 3D images 
% into a matrix of in-mask voxels. 
% The Coordinates of the in-mask voxels are returned as well as the
% paths of the input and mask images as they might be useful later.
% This step can be skipped if the data is already in a correct matrix
% format.
[Y, maskXYZ_vox, maskXYZ_mm, imagesPath, maskPath] = neuroCate_3Dimages2Matrix();
[nSubj, nVox] = size(Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Format/load the design matrices X and M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that no global intercept is automatically added. If needed, it must
% be added in X as is done below
load('covariatesForDemo')
X = [ones(nSubj,1), genderForDemo, adjGestationalAgeForDemo - mean(adjGestationalAgeForDemo)];

load('adjMethylationForDemo.mat') % contained M with 3 methylation loci
M = adjMethylationForDemo;
nML = size(M,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: estimate the number of unknown covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nZ =  zeros(1, nML);
nZMax = 10; % maximum number of unknown covariates allowed to be estimated
for iML = 1:nML
  nZ(iML) = neuroCate_estimateNumberOfFactors(Y, [X M(:,iML)], nZMax);
  % Note that additional information can be retrieved from estimateNumberOfFactors.
  % In the commented example below, the eigenvalues are returned in S, the
  % mock eigenvalue is in S0 (needed to allow nZ to be estimated to 0) and 
  % ER gives the Eigenvalue Ratio scores for nZ = 0, 1,...,nZMax. nZMax is
  % the maximum number of unknown covariates (the default value is 10).
  
  % [nZ(iML), S, S0, ER] = estimateNumberOfFactors(Y, [X M], nZMax);  
end
% select the most frequent estimate as a common nZ across all methylation loci
nZ_common = mode(nZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: compute a common resampling matrix across methylation locus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that this is not a mandatory step if nMN = 1. This step ensures that the
% same resampling matrix will be used for each methylation locus. This is
% useful if the jobs are split accross several computers/servers, which
% could generate different random numbers even with the same seed for the
% random generator
nB = 999; % this number can be reduced for the demo purpose if you are in a hurry
% set the seed of the random generator for repeatablity
seed = 0;
% compute a Wild Bootstrap resampling matrix  
resamplingMatrixWB = neuroCate_resamplingMatrix(nSubj, nB, 'wild bootstrap', seed);
% or compute a Permutation resampling matrix 
resamplingMatrixPerm = neuroCate_resamplingMatrix(nSubj, nB, 'permutation', seed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: run CATE at each methylation locus and build the max. global score
% distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxGlobalScoreF = zeros(nB + 1,1);
for iML = 1:nML
  % run CATE with Wild Bootstrap
  result = neuroCate_cate(Y, X, M(:,iML), nZ_common, 'inference', 'wb', 'nB', nB, 'resamplingMatrix', resamplingMatrixWB);
  % commented below, run CATE with permutation instead
  % results = neuroCate_cate(Y, X, M(:,iML), nZ, 'inference', 'perm', 'nB', nB, 'resamplingMatrix', resamplingMatrixPerm);
  % Additional options could be set in neuroCate_cate (for example instead 
  % of feeding a resampling matrix, a seed could be set instead to generate 
  % the resampling matrix inside the function; this option is more useful when
  % nML is 1)
  
  % save the results for each meth. locus on disk. 
  % Note that the variable "result" is a structure array containing a lot of fields of 
  % variables. If needed, save only the relevant fields to spare some disk memory
  str = sprintf('result_%i', iML); % please change the naming by something relevant for your study
  save(str, 'result');
  
  % update the maximum global statistic
  maxGlobalScoreF = max(maxGlobalScoreF, result.maxScoreF);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: compute the global FWER-corrected p-values for all original scores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minPGlobalFWERPerML = ones(1, nML); % minimum global FWER-corected p-value per meth. locus.
for iML = 1:nML
  str = sprintf('result_%i.mat', iML);  
  load(str);
  pGlobalFWER = ones(1, nVox);
  for b = 2:(nB+1)
    pGlobalFWER = pGlobalFWER + 1 .* (maxGlobalScoreF(b) >= result.scoreF);
  end
  pGlobalFWER = pGlobalFWER/ (nB  + 1);
  minPGlobalFWERPerML(iML) = min(pGlobalFWER);
  % add the global FWER-corected p-values to the results and save them
  % again - could be skipped if not useful 
  result.pGlobalFWER = pGlobalFWER;
  save(str, 'result');
  
  % produce output images here instead of in STEP 7 (see below for more 
  % details about this line of code)
  % This step can produce a lot of data and could be skipped or used only if, for example, 
  % min(pGlobalFWER) <= pThreshold or min(pGlobalFWER) < 1  
  neuroCate_createOutputImage(imagesPath(1,:),...
                              sprintf('lFWEp_neuroCate_FA_2Z_wb_ML%i.img', iML),...
                              -log10(pGlobalFWER),...
                              maskXYZ_vox,...
                              sprintf('-log10(FWER-corrected p-values) for ML #%i', iML))
end
% save the minimum global FWER-corected p-value per meth. locus. This is
% useful to quickly locate significant meth. loci. 
save('minPGlobalFWERPerML.mat', minPGlobalFWERPerML)

% If images info are saved, the results can be mapped into images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 7: map the outputs into images 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This step can also be done in the loop above if this is appropriate (See
% the line of code in Step 6). 
% We provide a function "neuroCate_createOutputImage" that uses SPM functions 
% for this step. This function map some data onto an image based on a
% template image (typically, one input image from the analysis).
% For p-values, we recommend to transform them into -log10(p-values) before 
% mapping them into images as it is easier to inspect visually than untransformed p-values.
% Below, we show only how to do it for -log10(p-values), but other results
% can be saved into images in a similar way.

% path of the template image 
templateImage = imagesPath(1,:);

% path of the output image (give the full path if the file need to be saved
% in another directory than the current working directory) 
outpuImage = sprintf('lFWEp_neuroCate_FA_2Z_wb_ML%i.img', iML);

% data for the output image 
load('result_1.mat');
data = -log10(result.pGlobalFWER);

% meaningful description of the output image
outputDescription = '-log10(FWER-corrected p-values) for ML #1';

% produce the output image
neuroCate_createOutputImage(templateImage, outpuImage, data, maskXYZ_vox, outputDescription)

