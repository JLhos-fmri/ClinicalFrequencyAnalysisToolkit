function CFA_RunAllJob(SetUpPara,Template,Para)
Indir = SetUpPara.Indir;
Outdir = SetUpPara.Outdir;
filp = which('CFA_RunAllJob.m');
[pth,nam,ext] = fileparts(filp);
spmjobdir = [pth,filesep,'NormJobmat',filesep];
Atlasdir = [pth,filesep,'TemplateUsed',filesep];
voxelsize = SetUpPara.OutVsize;
BrainMaskdir = [pth,filesep,'brainmask.nii'];

%% Normalize
if SetUpPara.NormOpt
    outdir = [Outdir,filesep,'Normalized',filesep];
    mkdir(outdir);
    sfold = dir(Indir);
    if SetUpPara.Normtype1 % epimode
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'image.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            outname = sfold(i+2).name;
            ClincFreq_Norm2MNI_EPI(indat,inroi,outdir,outname,voxelsize,spmjobdir)
        end
    elseif SetUpPara.Normtype2 % T1mode
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'image.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            outname = sfold(i+2).name;
            %             ClincFreq_Norm2MNI_T1(indat,inroi,outdir,outname,voxelsize,spmjobdir)
            ClincFreq_Norm2MNI_T1_clinc(indat,inroi,outdir,outname,voxelsize,spmjobdir)
        end
    elseif SetUpPara.Normtype3 % with T1mode
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'OtherMode.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            indatT1 = [Indir,filesep,sfold(i+2).name,filesep,'T1.nii'];
            outname = sfold(i+2).name;
            %             ClincFreq_Norm2MNI_withT1(indat,inroi,indatT1,outdir,outname,voxelsize,spmjobdir)
            ClincFreq_Norm2MNI_withT1_clinc(indat,inroi,indatT1,outdir,outname,voxelsize,spmjobdir)
        end
    elseif SetUpPara.Normtype4 % T2mode
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'image.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            outname = sfold(i+2).name;
            ClincFreq_Norm2MNI_T2(indat,inroi,outdir,outname,voxelsize,spmjobdir)
        end
    elseif SetUpPara.Normtype5 % CTmode
        parfor i = 1:length(sfold)-2
            indat = [Indir,filesep,sfold(i+2).name,filesep,'image.nii'];
            inroi = {[Indir,filesep,sfold(i+2).name,filesep,'ROI.nii']};
            outname = sfold(i+2).name;
            %             ClincFreq_Norm2MNI_CT(indat,inroi,outdir,outname,voxelsize,spmjobdir)
            ClincFreq_Norm2MNI_CT_clinc(indat,inroi,outdir,outname,voxelsize,spmjobdir)
        end
    end
    
    INDIR = [Outdir,filesep,'NormROI',filesep];
    mkdir(INDIR);
    copyfile([outdir,'ROI*.nii'],INDIR);
else
    INDIR = [Indir,filesep];
end

% dynamicBC_Reslice(AtlasTemp{i,1},AtlasTempU{i,1},voxelsize,0,[INDIR,filesep,ROIfil(1).name]);
Rfile = dir([INDIR,filesep,'*.nii']);
dynamicBC_Reslice(BrainMaskdir,[Outdir,filesep,'Bmask.nii'],voxelsize*[1 1 1],0,[INDIR,filesep,Rfile(1).name]);
[vbmask dbmask0] = Dynamic_read_dir_NIFTI([Outdir,filesep,'Bmask.nii']);
dbmask = reshape(dbmask0,vbmask.dim);
%%
atlasnum = 0;
if Template.Fivetissue
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_FiveTissue.nii'];
    load([Atlasdir,'AtlasU_FiveTissue.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'FiveTissue';
end
if Template.FivetissueLR
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_FiveTissueLR.nii'];
    load([Atlasdir,'AtlasU_FiveTissueLR.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'FiveTissueLR';
end
if Template.JHUtotal
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_JHUtotal.nii'];
    load([Atlasdir,'AtlasU_JHUtotal.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'JHUtotal';
end
if Template.LPBA40
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_LPBA40.nii'];
    load([Atlasdir,'AtlasU_LPBA40.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'LPBA40';
end
if Template.AAL
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_AAL.nii'];
    load([Atlasdir,'AtlasU_AAL.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'AAL';
end
if Template.HOAcortex
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_HOAcortex.nii'];
    load([Atlasdir,'AtlasU_HOAcortex.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'HOAcortex';
end
if Template.HOAsubcortex
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_HOAsubcortex.nii'];
    load([Atlasdir,'AtlasU_HOAsubcortex.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'HOAsubcortex';
end
if Template.MesulamHOA
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_MesulamHOA.nii'];
    load([Atlasdir,'AtlasU_MesulamHOA.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'MesulamHOA';
end
if Template.JHUwhite
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_JHUwhite.nii'];
    load([Atlasdir,'AtlasU_JHUwhite.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'JHUwhite';
end
if Template.ICBMwhite
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_ICBMwhite.nii'];
    load([Atlasdir,'AtlasU_ICBMwhite.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'ICBMwhite';
end
if Template.Ventricle
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = [Atlasdir,'AtlasU_Ventricle.nii'];
    load([Atlasdir,'AtlasU_Ventricle.mat']);
    AtlasTemp{atlasnum,2} = ROItab;
    AtlasTemp{atlasnum,3} = 'Ventricle';
end
if Template.Other
    atlasnum = atlasnum+1;
    AtlasTemp{atlasnum,1} = Template.OtherROI;
    load(Template.OtherROIname);
    AtlasTemp{atlasnum,2} = ROItab;
    [pths,nams,exts] = fileparts(Template.OtherROI);
    AtlasTemp{atlasnum,3} = nams;
end

%%
if ~isempty(whos('AtlasTemp'))
    OutdirAtlas = [Outdir,filesep,'ReslicedAtlas',filesep];
    mkdir(OutdirAtlas)
    ROIfil = dir([INDIR,filesep,'*.nii']);
    for i = 1:size(AtlasTemp,1)
        tempdir = AtlasTemp{i,1};
        [pt,na,ex] = fileparts(tempdir);
        if Template.Other&&i==size(AtlasTemp,1)
            AtlasTempU{i,1} = [OutdirAtlas,filesep,'AtlasU_',na,ex];
        else            
            AtlasTempU{i,1} = [OutdirAtlas,filesep,na,ex];
        end
        AtlasTempU{i,2} = AtlasTemp{i,2};
        AtlasTempU{i,3} = AtlasTemp{i,3};
        dynamicBC_Reslice(AtlasTemp{i,1},AtlasTempU{i,1},voxelsize,0,[INDIR,filesep,ROIfil(1).name]);
    end
    
    for i = 1:size(AtlasTempU,1)
        [VT DT] = Dynamic_read_dir_NIFTI(AtlasTempU{i,1});
        ATMark = AtlasTempU{i,2};
        for j = 1:size(ATMark,1)
            labval(j,1) = ATMark{j,1};
            labname{j,1} = ATMark{j,2};
        end
        indTab = find(labval>0);
        AtlasTempU{i,4} = labval(indTab);
        AtlasTempU{i,5} = labname(indTab);
        for j = 1:length(indTab)
            AtlasTempInd{j,1} = find(DT==labval(indTab(j)));
        end
        AtlasTempU{i,6} = AtlasTempInd;
        clear labval labname AtlasTempInd
    end
    if atlasnum==0
        AtlasTempU = [];
    end
    
    if Para.Heat.V||Para.Heat.R
        OutHeatdir = [Outdir,filesep,'HeatNumber',filesep];
        mkdir(OutHeatdir);
        disp('Now Calculate HeatNumber!')
        ClincFreq_CalHeatNum(INDIR,OutHeatdir,SetUpPara,Para,AtlasTempU,dbmask)
        % stat needed
        ClincFreq_HeatNumStat(INDIR,OutHeatdir,SetUpPara,Para,AtlasTempU);
        disp('HeatNumber Calculation finished')
    end
    if Para.CentNum.R||Para.WeiCent.V||Para.WeiCent.R
        outdirCent = [Outdir,filesep,'ForCentRelated',filesep];
        if isempty(dir(outdirCent))
            mkdir(outdirCent)
            ClincFreq_CentPointForVolume_pre(INDIR,outdirCent,Para,spmjobdir);
        else
            if Para.WeiCent.V||Para.WeiCent.R
                if isempty(dir([outdirCent,'VolWeiCent']))
                    ClincFreq_CentPointForVolume_pre(INDIR,outdirCent,Para,spmjobdir);
                end
            end
        end
    end
    if Para.CentNum.R
        OutCentdir = [Outdir,filesep,'CenterNumber',filesep];
        mkdir(OutCentdir);
        disp('Now Calculate CenterNumber!')
        ClincFreq_CalCentNum(outdirCent,OutCentdir,SetUpPara,Para,AtlasTempU)
        % stat needed
        ClincFreq_CentNumStat(outdirCent,OutCentdir,SetUpPara,Para,AtlasTempU);
        disp('CenterNumber finished')
    end
    if Para.AffVolume.R
        OutAffVoldir = [Outdir,filesep,'AffectedVolume',filesep];
        mkdir(OutAffVoldir);
        disp('Now Calculate AffectedVolume!')
        ClincFreq_AffVolume(INDIR,OutAffVoldir,SetUpPara,Para,AtlasTempU)
        % stat needed
        ClincFreq_AffVolStat(INDIR,OutAffVoldir,SetUpPara,Para,AtlasTempU);
        disp('AffectedVolume finished')
    end
    if Para.WeiCent.V||Para.WeiCent.R
        OutWeiCentdir = [Outdir,filesep,'VolumeWeightedCenter',filesep];
        mkdir(OutWeiCentdir);
        disp('Now Calculate VolumeWeightedCenter!')
        ClincFreq_CalWeiCent(outdirCent,OutWeiCentdir,SetUpPara,Para,AtlasTempU,dbmask)
        % stat needed
        ClincFreq_WeiCentStat(outdirCent,OutWeiCentdir,SetUpPara,Para,AtlasTempU);
        disp('VolumeWeightedCenter finished')
    end
end

end