function ClincFreq_Norm2MNI_T1_clinc(indat,inroi,outdir,outname,voxelsize,pth)
[patdat,namdat,extdat] = fileparts(indat);
% [patroi,namroi,extroi] = fileparts(inroi);
% load([pth,filesep,'Norm2T1_segment.mat']);
load([pth,filesep,'Norm2T1_clincal.mat']);

pathspm = which('spm.m');
[patspm,namspm,extspm] = fileparts(pathspm);

MATLABBAT_seg = matlabbatch;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.anat{1} = indat;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.les{1} = inroi{1,1};
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.vox = [1,1,1]*voxelsize;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.bb = [-90 -126 -72;90 90 108];
spm_jobman('run',MATLABBAT_seg)

copyfile([patdat,filesep,'w',namdat,extdat],[outdir,filesep,'image_',outname,'_w',num2str(voxelsize),'_',namdat,extdat])

[patroi,namroi,extroi] = fileparts(inroi{1,1});
copyfile([patroi,filesep,'ws',namroi,extroi],[outdir,filesep,'ROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi])

end