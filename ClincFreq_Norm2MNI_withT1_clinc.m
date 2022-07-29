function ClincFreq_Norm2MNI_withT1_clinc(indat,inroi,indatT1,outdir,outname,voxelsize,pth)
[patdat,namdat,extdat] = fileparts(indat);
% [patroi,namroi,extroi] = fileparts(inroi);
[patdatT1,namdatT1,extdatT1] = fileparts(indatT1);
%%

[patdat,namdat,extdat] = fileparts(indat);
% [patroi,namroi,extroi] = fileparts(inroi);
% load([pth,filesep,'Norm2T1_segment.mat']);
load([pth,filesep,'Norm2T1_clincal.mat']);

pathspm = which('spm.m');
[patspm,namspm,extspm] = fileparts(pathspm);

MATLABBAT_seg = matlabbatch;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.anat{1} = indatT1;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.les{1} = inroi{1,1};
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.t2{1} = indat;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.vox = [1,1,1]*voxelsize;
MATLABBAT_seg{1}.spm.tools.MRI.MRnormseg.bb = [-90 -126 -72;90 90 108];
spm_jobman('run',MATLABBAT_seg)

% copyfile([patdat,filesep,'w',namdat,extdat],[outdir,filesep,'image_',outname,'_w',num2str(voxelsize),'_',namdat,extdat])
copyfile([patdat,filesep,'w',namdatT1,extdatT1],[outdir,filesep,'image_T1_',outname,'_w',num2str(voxelsize),'_',namdat,extdat])

[patroi,namroi,extroi] = fileparts(inroi{1,1});
copyfile([patroi,filesep,'wsr',namroi,extroi],[outdir,filesep,'ROI_',outname,'_w',num2str(voxelsize),'_',namroi,extroi])


end
