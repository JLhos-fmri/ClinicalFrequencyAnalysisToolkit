function CFA_ReOrient_Twomode(T1dir,Otherdir,ROIdir)
%%
coivox = ones(4,1); %center of intensity

% hdr = spm_vol([Imagedir,',1']); %load header
% img = spm_read_vols(hdr); %load image data
[hdr,imgr] = Dynamic_read_dir_NIFTI(T1dir);
img = reshape(imgr,hdr.dim);
img = img - min(img(:));
img(isnan(img)) = 0;
sumTotal = sum(img(:));
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZ_mm = hdr.mat * coivox; %convert from voxels to millimeters
hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);

hdr.private.mat = hdr.mat;
hdr.private.mat0 = hdr.mat;

DynamicBC_write_NIFTI(img,hdr,T1dir);
%%

coivox = ones(4,1); %center of intensity

% hdr = spm_vol([Imagedir,',1']); %load header
% img = spm_read_vols(hdr); %load image data
[hdr,imgr] = Dynamic_read_dir_NIFTI(Otherdir);
img = reshape(imgr,hdr.dim);
img = img - min(img(:));
img(isnan(img)) = 0;
sumTotal = sum(img(:));
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZ_mm = hdr.mat * coivox; %convert from voxels to millimeters
hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);

hdr.private.mat = hdr.mat;
hdr.private.mat0 = hdr.mat;

DynamicBC_write_NIFTI(img,hdr,Otherdir);

[hdrROI,imgROI] = Dynamic_read_dir_NIFTI(ROIdir);
imgROI2 = reshape(imgROI,hdrROI.dim);
imgROI2(isnan(imgROI2)) = 0;
imgROI2(imgROI2>0) = 1;
DynamicBC_write_NIFTI(imgROI2,hdr,ROIdir);


end