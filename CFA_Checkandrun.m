function CFA_Checkandrun(H)
% function CFA_Checkandrun
pathss = which('CFA_Checkandrun.m');
[pth nam ext] = fileparts(pathss);
IOsetupdir = [pth,filesep,'IOparameter',filesep];
iosetup = dir('SetUpPara.mat');
if isempty(iosetup)
    iosetup = dir([IOsetupdir,'SetUpPara.mat']);
    if isempty(iosetup)
        [fil pg] = uigetfile('*.mat','IO mats (SetUpPara.mat)');
        load([pg,filesep,fil]);
    else
        load([IOsetupdir,'SetUpPara.mat']);
    end
else
    load('SetUpPara.mat');
end
TempAndPara = dir('Setup_TempAndPara.mat');
if isempty(iosetup)
    iosetup = dir([IOsetupdir,'Setup_TempAndPara.mat']);
    if isempty(iosetup)
        [fil pg] = uigetfile('*.mat','Template&Parameter mats (Setup_TempAndPara.mat)');
        load([pg,filesep,fil]);
    else
        load([IOsetupdir,'Setup_TempAndPara.mat']);
    end
else
    load('Setup_TempAndPara.mat');
end
%%
HCheckRun.fig = figure('unit','norm',...
    'pos',[0.25,0.2,0.5,0.6],...
    'menubar','none',...
    'color','w',...
    'name','Check and Run');

HCheckRun.Reset = uicontrol('parent',HCheckRun.fig,...
    'unit','norm',...
    'pos',[0.2,0.1,0.2,0.1],...
    'style','pushbutton',...
    'string','Reset');
HCheckRun.Run = uicontrol('parent',HCheckRun.fig,...
    'unit','norm',...
    'pos',[0.6,0.1,0.2,0.1],...
    'style','pushbutton',...
    'string','Run');

HCheckRun.NormOpt = uibuttongroup('parent',HCheckRun.fig,...
    'unit','norm',...
    'pos',[0.1,0.6,0.8,0.35],...
    'title','I/O and Normalize Operation');
HCheckRun.TempOpt = uibuttongroup('parent',HCheckRun.fig,...
    'unit','norm',...
    'pos',[0.1,0.45,0.8,0.15],...
    'title','Template Selection');
HCheckRun.ParaOpt = uibuttongroup('parent',HCheckRun.fig,...
    'unit','norm',...
    'pos',[0.1,0.25,0.8,0.2],...
    'title','Parameter Selection');

%%
SHOWINFO.IN = uicontrol('parent',HCheckRun.NormOpt,...
    'unit','norm',...
    'pos',[0.01,0.91,0.98,0.08],...
    'style','text',...
    'string',['Indir: ',SetUpPara.Indir]);
% SetUpPara.Indir
SHOWINFO.OUT = uicontrol('parent',HCheckRun.NormOpt,...
    'unit','norm',...
    'pos',[0.01,0.81,0.98,0.08],...
    'style','text',...
    'string',['Outdir:',SetUpPara.Outdir]);
% SetUpPara.Outdir
if SetUpPara.NormOpt
    if SetUpPara.Normtype1
        SHOWINFO.Norm = uicontrol('parent',HCheckRun.NormOpt,...
            'unit','norm',...
            'pos',[0.01,0.71,0.98,0.08],...
            'style','text',...
            'string',['Need Normalize: EPI mode, VoxelSize: ',num2str(SetUpPara.OutVsize)]);
    elseif SetUpPara.Normtype2
        SHOWINFO.Norm = uicontrol('parent',HCheckRun.NormOpt,...
            'unit','norm',...
            'pos',[0.01,0.71,0.98,0.08],...
            'style','text',...
            'string',['Need Normalize: T1 mode, VoxelSize: ',num2str(SetUpPara.OutVsize)]);
    elseif SetUpPara.Normtype3
        SHOWINFO.Norm = uicontrol('parent',HCheckRun.NormOpt,...
            'unit','norm',...
            'pos',[0.01,0.71,0.98,0.08],...
            'style','text',...
            'string',['Need Normalize: with T1 mode, VoxelSize: ',num2str(SetUpPara.OutVsize)]);
    elseif SetUpPara.Normtype4
        SHOWINFO.Norm = uicontrol('parent',HCheckRun.NormOpt,...
            'unit','norm',...
            'pos',[0.01,0.71,0.98,0.08],...
            'style','text',...
            'string',['Need Normalize: T2 mode, VoxelSize: ',num2str(SetUpPara.OutVsize)]);
    elseif SetUpPara.Normtype5
        SHOWINFO.Norm = uicontrol('parent',HCheckRun.NormOpt,...
            'unit','norm',...
            'pos',[0.01,0.71,0.98,0.08],...
            'style','text',...
            'string',['Need Normalize: CT mode, VoxelSize: ',num2str(SetUpPara.OutVsize)]);
    end
else    
    SHOWINFO.Norm = uicontrol('parent',HCheckRun.NormOpt,...
        'unit','norm',...
        'pos',[0.01,0.71,0.98,0.08],...
        'style','text',...
        'string','No Normalize Operation');
end
% SetUpPara.NormOpt
HCheckRun.Checklist = uicontrol('parent',HCheckRun.NormOpt,...
    'unit','norm',...
    'pos',[0.01,0.01,0.28,0.68],...
    'style','list');

% SetUpPara.ParaOut.LabName
% SetUpPara.ParaOut.datalab
% SetUpPara.ParaOut.Excludedlab
% SetUpPara.ParaOut.Labtxt
for i = 1:length(SetUpPara.ParaOut)
    listname{i,1} = SetUpPara.ParaOut(i).LabName;
end
set(HCheckRun.Checklist,'string',listname);
set(HCheckRun.Checklist,'val',1);

HCheckRun.LabInfo = uibuttongroup('parent',HCheckRun.NormOpt,...
    'unit','norm',...
    'pos',[0.31,0.01,0.68,0.68],...
    'title','Label Information');
HCheckRun.Labtxt = uicontrol('parent',HCheckRun.LabInfo,...
    'unit','norm',...
    'pos',[0.01,0.76,0.98,0.23],...
    'style','text',...
    'string',['Label from: ',SetUpPara.ParaOut(1).Labtxt]);
HCheckRun.Labthr = uicontrol('parent',HCheckRun.LabInfo,...
    'unit','norm',...
    'pos',[0.01,0.51,0.98,0.23],...
    'style','text',...
    'string',['Cutoff Threshold: ', mat2str(SetUpPara.ParaOut(1).thr)]);
HCheckRun.Label = uicontrol('parent',HCheckRun.LabInfo,...
    'unit','norm',...
    'pos',[0.01,0.26,0.98,0.23],...
    'style','text');
labStrings = [num2str(length(SetUpPara.ParaOut(1).datalab)),' periods exist:'];
for i = 1:length(SetUpPara.ParaOut(1).datalab)
    labStrings = [labStrings,' ',num2str(length(SetUpPara.ParaOut(1).datalab{i}))];
end
set(HCheckRun.Label,'String',labStrings)
HCheckRun.Excluded = uicontrol('parent',HCheckRun.LabInfo,...
    'unit','norm',...
    'pos',[0.01,0.01,0.98,0.23],...
    'style','text',...
    'string',['Excluded Number: ',num2str(length(SetUpPara.ParaOut(1).Excludedlab))]);
%%
% Template HCheckRun.TempOpt
SHOWINFO.Temp1 = uicontrol('parent',HCheckRun.TempOpt,...
    'unit','norm',...
    'pos',[0.01,0.51,0.98,0.48],...
    'style','text');
SHOWINFO.Temp2 = uicontrol('parent',HCheckRun.TempOpt,...
    'unit','norm',...
    'pos',[0.01,0.01,0.98,0.48],...
    'style','text');
StringTemp1 = 'Atlas:';
if Template.Fivetissue
    StringTemp1 = [StringTemp1,' Five Tissue;'];
end
if Template.FivetissueLR
    StringTemp1 = [StringTemp1,' Five Tissue(LR);'];
end
if Template.JHUtotal
    StringTemp1 = [StringTemp1,' JHU(Gray+White);'];
end
if Template.LPBA40
    StringTemp1 = [StringTemp1,' LPBA40;'];
end
if Template.AAL
    StringTemp1 = [StringTemp1,' AAL;'];
end
if Template.HOAcortex
    StringTemp1 = [StringTemp1,' HOA(cortex);'];
end
if Template.HOAsubcortex
    StringTemp1 = [StringTemp1,' HOA(subcortex);'];
end
if Template.MesulamHOA
    StringTemp1 = [StringTemp1,' HOA(Mesulam);'];
end
if Template.JHUwhite
    StringTemp1 = [StringTemp1,' JHU(white);'];
end
if Template.ICBMwhite
    StringTemp1 = [StringTemp1,' ICBM(white);'];
end
if Template.Ventricle
    StringTemp1 = [StringTemp1,' Ventricle;'];
end
set(SHOWINFO.Temp1,'string',StringTemp1);
if Template.Other
    StringTemp2 = ['Private ROI: ', Template.OtherROIname,'; with NIfTI: ',Template.OtherROI];
else
    StringTemp2 = 'No private ROI selection';
end
set(SHOWINFO.Temp2,'string',StringTemp2);
% %%  HCheckRun.ParaOpt
% Para.Heat
% Para.CentNum
% Para.AffVolume
% Para.WeiCent
SHOWINFO.Heat = uibuttongroup('parent',HCheckRun.ParaOpt,...
    'unit','norm',...
    'pos',[0/5,0,1/5,1],...
    'title','Heat Number');
% Para.Heat
% Heatstring
heatnum = 0;
try
    Para = ParaSet;
catch
    
end
if Para.Heat.V
    heatnum = heatnum+1;
    Heatstring{heatnum,1} = 'Voxel-wise analysis';
end
if Para.Heat.R
    heatnum = heatnum+1;
    Heatstring{heatnum,1} = ['ROI-wise analysis: ROI > ',num2str(Para.Heat.Rper),'%'];
end
if Para.Heat.StChisquare
    heatnum = heatnum+1;
    Heatstring{heatnum,1} = 'Use Chi-square Test';
end
if Para.Heat.StPermutation
    heatnum = heatnum+1;
    Heatstring{heatnum,1} = 'Use Permutation Test';
end
for i = 1:heatnum
    Showheattemp(i) = uicontrol('parent',SHOWINFO.Heat,...
        'unit','norm',...
        'pos',[0,1-(i/4),1,1/4],...
        'style','text',...
        'string',Heatstring{i,1});    
end

SHOWINFO.Cent = uibuttongroup('parent',HCheckRun.ParaOpt,...
    'unit','norm',...
    'pos',[1/5,0,1/5,1],...
    'title','Center Number');
% Para.CentNum

centnum = 0;
if Para.CentNum.R
    centnum = centnum+1;
    Centstring{centnum,1} = 'ROI-wise analysis';
end
if Para.Heat.StChisquare
    centnum = centnum+1;
    Centstring{centnum,1} = 'Use Chi-square Test';
end
if Para.Heat.StPermutation
    centnum = centnum+1;
    Centstring{centnum,1} = 'Use Permutation Test';
end
for i = 1:centnum
    Showheattemp(i) = uicontrol('parent',SHOWINFO.Cent,...
        'unit','norm',...
        'pos',[0,1-(i/4),1,1/4],...
        'style','text',...
        'string',Centstring{i,1});    
end

SHOWINFO.AffVol = uibuttongroup('parent',HCheckRun.ParaOpt,...
    'unit','norm',...
    'pos',[2/5,0,1/5,1],...
    'title','Affected Volume');
% Para.AffVolume

Affvol = 0;
if Para.AffVolume.R
    Affvol = Affvol+1;
    Affvolstring{Affvol,1} = 'ROI-wise analysis';
end
if Para.AffVolume.StPermutation
    Affvol = Affvol+1;
    Affvolstring{Affvol,1} = 'Use Permutation Test';
end
for i = 1:Affvol
    Showheattemp(i) = uicontrol('parent',SHOWINFO.AffVol,...
        'unit','norm',...
        'pos',[0,1-(i/4),1,1/4],...
        'style','text',...
        'string',Affvolstring{i,1});    
end

SHOWINFO.WeiCent = uibuttongroup('parent',HCheckRun.ParaOpt,...
    'unit','norm',...
    'pos',[3/5,0,1/5,1],...
    'title','Weighted Centre');
% Para.Cent

WeiCentNum = 0;
if Para.WeiCent.V
    WeiCentNum = WeiCentNum+1;
    WeiCentstring{WeiCentNum,1} = 'Voxel-wise analysis';
end
if Para.WeiCent.R
    WeiCentNum = WeiCentNum+1;
    WeiCentstring{WeiCentNum,1} = 'ROI-wise analysis';
end
if Para.WeiCent.StPermutation
    WeiCentNum = WeiCentNum+1;
    WeiCentstring{WeiCentNum,1} = 'Use Permutation Test';
end
for i = 1:WeiCentNum
    Showheattemp(i) = uicontrol('parent',SHOWINFO.WeiCent,...
        'unit','norm',...
        'pos',[0,1-(i/4),1,1/4],...
        'style','text',...
        'string',WeiCentstring{i,1});    
end
%%
set(HCheckRun.Checklist,'callback',{@CheckListLab,HCheckRun,SetUpPara})
set(HCheckRun.Reset,'callback',{@ResetFun,HCheckRun,H});
set(HCheckRun.Run,'callback',{@RunJob,HCheckRun,H,SetUpPara,Template,Para})
end

function CheckListLab(varargin)
HCheckRun = varargin{3};
SetUpPara = varargin{4};
Vals = get(HCheckRun.Checklist,'val');

set(HCheckRun.Labtxt,'string',['Label from: ',SetUpPara.ParaOut(Vals).Labtxt]);
set(HCheckRun.Labthr,'string',['Cutoff Threshold: ', mat2str(SetUpPara.ParaOut(Vals).thr)]);

labStrings = [num2str(length(SetUpPara.ParaOut(Vals).datalab)),' periods exist:'];
for i = 1:length(SetUpPara.ParaOut(Vals).datalab)
    labStrings = [labStrings,' ',num2str(length(SetUpPara.ParaOut(Vals).datalab{i}))];
end
set(HCheckRun.Label,'String',labStrings)

set(HCheckRun.Excluded,'string',['Excluded Number: ',num2str(length(SetUpPara.ParaOut(Vals).Excludedlab))]);
end
function ResetFun(varargin)
HCheckRun = varargin{3};
H = varargin{4};
close(HCheckRun.fig);
close(H.fig);
% CFATmain;
% CFATmainVer1;
CFAT;
end
function RunJob(varargin)
HCheckRun = varargin{3};
H = varargin{4};
SetUpPara = varargin{5};
Template = varargin{6};
Para = varargin{7};

CFA_RunAllJob(SetUpPara,Template,Para);
uiwait(msgbox('Finished!'));
close(HCheckRun.fig);
end