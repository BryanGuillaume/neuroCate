function  neuroCate_createOutputImage(templateImage, outpuImage, data, maskXYZ_vox, outputDescription)

VO = spm_vol(templateImage);
VO.fname = outpuImage;
VO.descrip = outputDescription;

% create the output image
VO = spm_create_vol(VO);

% create a 3D array mask
Q = cumprod([1,VO.dim(1:2)])* maskXYZ_vox - ...
  sum(cumprod(VO.dim(1:2)));

% create a 3D array of nan
tmp = nan(VO.dim);

% fill in the in-mask voxels
tmp(Q) = data;

% write the output image
spm_write_vol(VO, tmp);
