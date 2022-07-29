function ClincFreq_Norm2MNI_T2(indat,inroi,outdir,outname,voxelsize,pth)
[patdat,namdat,extdat] = fileparts(indat);
% [patroi,namroi,extroi] = fileparts(inroi);
load([pth,filesep,'Norm2EPI_est.mat']);
pathspm = which('spm.m');
[patspm,namspm,extspm] = fileparts(pathspm);

MATLABBAT_est = matlabbatch;
MATLABBAT_est{1}.spm.tools.oldnorm.est.subj.source{1} = indat;
MATLABBAT_est{1}.spm.tools.oldnorm.est.eoptions.template{1} = fullfile(patspm,'toolbox','OldNorm','T2.nii');

spm_jobman('run',MATLABBAT_est);

load([pth,filesep,'Norm2EPIandCT_writeImage.mat'])
MATLABBAT_writeimage = matlabbatch;
MATLABBAT_writeimage{1}.spm.tools.oldnorm.write.subj.matname{1} = fullfile(patdat,[namdat,'_sn.mat']);
MATLABBAT_writeimage{1}.spm.tools.oldnorm.write.subj.resample{1} = indat;
MATLABBAT_writeimage{1}.spm.tools.oldnorm.write.roptions.vox = [1,1,1]*voxelsize;
MATLABBAT_writeimage{1}.spm.tools.oldnorm.write.roptions.prefix = ['w',num2str(voxelsize)];
spm_jobman('run',MATLABBAT_writeimage);

load([pth,filesep,'Norm2EPIandCT_writeROI.mat'])
MATLABBAT_writeroi = matlabbatch;
MATLABBAT_writeroi{1}.spm.tools.oldnorm.write.subj.matname{1} = fullfile(patdat,[namdat,'_sn.mat']);
for i = 1:length(inroi)
    MATLABBAT_writeroi{1}.spm.tools.oldnorm.write.subj.resample{i,1} = inroi{i,1};
end
MATLABBAT_writeroi{1}.spm.tools.oldnorm.write.roptions.vox = [1,1,1]*voxelsize;
MATLABBAT_writeroi{1}.spm.tools.oldnorm.write.roptions.prefix = ['w',num2str(voxelsize)];
spm_jobman('run',MATLABBAT_writeroi);

copyfile([patdat,filesep,'w',num2str(voxelsize),namdat,extdat],[outdir,filesep,'image_',outname,'_w',num2str(voxelsize),'_',namdat,extdat]);
% copyfile([patroi,filesep,'w',num2str(voxelsize),namroi,extdat],[outdir,filesep,'ROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi]);
for i = 1:length(inroi)
    [patroi,namroi,extroi] = fileparts(inroi{i,1});
    copyfile([patroi,filesep,'w',num2str(voxelsize),namroi,extdat],[outdir,filesep,'ROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi]);
end
end