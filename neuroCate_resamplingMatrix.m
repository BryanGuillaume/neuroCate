function [resamplingMatrix] = neuroCate_resamplingMatrix(nSubj, nB, type, seed)

resamplingMatrix = zeros(nSubj, nB);
% set the seed
rng(seed);

if strcmpi(type, 'permutation') || strcmpi(type, 'perm')
  for iB = 1:nB
    resamplingMatrix(:,iB) = randperm(nSubj);
  end
elseif strcmpi(type, 'wild bootstrap') || strcmpi(type, 'wb')
      resamplingMatrix = randi([0 1], nSubj, nB)*2 - 1;
else
  error('unknown type of resampling matrix');
end
  

