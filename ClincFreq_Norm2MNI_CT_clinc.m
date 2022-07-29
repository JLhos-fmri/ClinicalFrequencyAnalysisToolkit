function ClincFreq_Norm2MNI_CT_clinc(indat,inroi,outdir,outname,voxelsize,pth)
%% Need to rebuild

[patdat,namdat,extdat] = fileparts(indat);
% [patroi,namroi,extroi] = fileparts(inroi);
load([pth,filesep,'norm2ct.mat']);
% load([pth,filesep,'Norm2CT_est.mat']);
pathspm = which('spm.m');
[patspm,namspm,extspm] = fileparts(pathspm);

% indat0 = indat;
% indat = [patdat,filesep,'c',namdat,extdat];
[patdat,namdat,extdat] = fileparts(indat);

MATLABBAT_est = matlabbatch;
MATLABBAT_est{1}.spm.tools.MRI.CTnorm.images{1} = indat;
% MATLABBAT_est{1}.spm.tools.MRI.CTnorm.bb
MATLABBAT_est{1}.spm.tools.MRI.CTnorm.vox = [1 1 1]*voxelsize;
MATLABBAT_est{1}.spm.tools.MRI.CTnorm.ctles{1} = inroi{1,1};
spm_jobman('run',MATLABBAT_est);

copyfile([patdat,filesep,'w',namdat,extdat],[outdir,filesep,'image_',outname,'_w',num2str(voxelsize),'_',namdat,extdat]);
[patroi,namroi,extroi] = fileparts(inroi{1,1});
copyfile([patroi,filesep,'ws',namroi,extdat],[outdir,filesep,'ROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi]);

end