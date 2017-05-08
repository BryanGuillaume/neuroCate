function [Y, maskXYZ_vox, maskXYZ_mm, imagesPath] = neuroCate_3Dimages2Matrix()

imagesPath = spm_select([1 Inf], 'image', 'Select 3D images to analyse');
maskPath = spm_select(1, 'image', 'Select a mask');
% maskPath = spm_select([0 1], 'image', 'Select a mask if needed (recommended)');

% if ~isempty(maskPath)
  VM = spm_vol(maskPath);
  [mask, maskXYZ_mm] = spm_read_vols(VM);
  maskXYZ_vox = round(VM.mat \ [maskXYZ_mm; ones(1,prod(VM.dim))]);
  maskXYZ_vox = maskXYZ_vox(1:3,:);
  maskXYZ_vox = maskXYZ_vox(:,mask(:) > 0);
  Y = spm_get_data(imagesPath, maskXYZ_vox);
% else
%   VI = spm_vol(imagesPath);
%  	Y = spm_read_vols(VI);
% end


