function [Y, maskXYZ_vox, maskXYZ_mm] = neuroCate_3Dimages2Matrix()

imagesPath = spm_select([1 Inf], 'image', 'Select 3D images to analyse');
maskPath = spm_select(1, 'image', 'Select a mask');
% maskPath = spm_select([0 1], 'image', 'Select a mask if needed (recommended)');

% if ~isempty(maskPath)
  VI = spm_vol(maskPath);
  [mask, maskXYZ_mm] = spm_read_vols(VI);
  maskXYZ_vox = round(VI.mat \ [maskXYZ_mm; ones(1,prod(VI.dim))]);
  maskXYZ_vox = maskXYZ_vox(1:3,:);
  maskXYZ_vox = maskXYZ_vox(:,mask(:) > 0);
  Y = spm_get_data(imagesPath, maskXYZ_vox);
% else
%   VI = spm_vol(imagesPath);
%  	Y = spm_read_vols(VI);
% end


