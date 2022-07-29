function CFA_inout(H)
% function CFA_inout
movegui(H.fig,'west')

pathss = which('CFA_inout.m');
[pth nam ext] = fileparts(pathss);
IOsetupdir = [pth,filesep,'IOparameter',filesep];
if ~isempty(dir(IOsetupdir))
    rmdir(IOsetupdir,'s')
    mkdir(IOsetupdir);
else
    mkdir(IOsetupdir);
end
load([pth,filesep,'colorforlab.mat']);
HIO.colorforlab = colorforlab;
HIO.IOsetupdir = IOsetupdir;
HIO.fig = figure('unit','norm',...
    'pos',[0.2,0.2,0.6,0.6],...
    'menubar','none',...
    'color','w',...
    'name','I/O & label setup');
HIO.Normsel = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.1,0.9,0.1,0.05],...
    'style','rad',...
    'string','Need Normalize');
HIO.withNorm = uibuttongroup('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.2,0.9,0.7,0.05],...
    'title','Normalize Options');
HIO.Normtype1 = uicontrol('parent',HIO.withNorm,...
    'unit','norm',...
    'pos',[0.01,0.01,0.18,0.98],...
    'style','rad',...
    'string','EPI mode');
HIO.Normtype2 = uicontrol('parent',HIO.withNorm,...
    'unit','norm',...
    'pos',[0.21,0.01,0.18,0.98],...
    'style','rad',...
    'string','T1 mode');
HIO.Normtype3 = uicontrol('parent',HIO.withNorm,...
    'unit','norm',...
    'pos',[0.41,0.01,0.18,0.98],...
    'style','rad',...
    'string','withT1 mode');
HIO.Normtype4 = uicontrol('parent',HIO.withNorm,...
    'unit','norm',...
    'pos',[0.61,0.01,0.18,0.98],...
    'style','rad',...
    'string','T2 mode');
HIO.Normtype5 = uicontrol('parent',HIO.withNorm,...
    'unit','norm',...
    'pos',[0.81,0.01,0.18,0.98],...
    'style','rad',...
    'string','CT mode');

%%
HIO.intext = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.1,0.85,0.05,0.04],...
    'style','text',...
    'string','Input Directory');
HIO.inedit = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.15,0.85,0.7,0.04],...
    'style','edit',...
    'string','NuLL');
HIO.inpb = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.85,0.85,0.05,0.04],...
    'style','pushbutton',...
    'string','...');
HIO.outtext = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.1,0.8,0.05,0.04],...
    'style','text',...
    'string','Output Directory');
HIO.outedit = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.15,0.8,0.7,0.04],...
    'style','edit',...
    'string','NuLL');
HIO.outpb = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.85,0.8,0.05,0.04],...
    'style','pushbutton',...
    'string','...');
%%
HIO.lablisttitle = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.1,0.725,0.25,0.05],...
    'style','text',...
    'string','Label Lists');


HIO.lablist = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.1,0.2,0.25,0.5],...
    'style','list');
%     'string',{'test';'test2';'test3'},...
%%
HIO.labOpt = uibuttongroup('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.4,0.2,0.2,0.575],...
    'title','Label Operation');

HIO.labOpt_add = uicontrol('parent',HIO.labOpt,...
    'unit','norm',...
    'pos',[0.05,0.9,0.425,0.05],...
    'style','pushbutton',...
    'String','Add Current List');
HIO.labOpt_del = uicontrol('parent',HIO.labOpt,...
    'unit','norm',...
    'pos',[0.525,0.9,0.425,0.05],...
    'style','pushbutton',...
    'String','Delete Current List');

HIO.labOpt_excluded = uibuttongroup('parent',HIO.labOpt,...
    'unit','norm',...
    'pos',[0.05,0.6,0.9,0.3],...
    'title','Exclude Condition');
HIO.labOpt_excluded_lessOpt = uicontrol('parent',HIO.labOpt_excluded,...
    'unit','norm',...
    'pos',[0.1,0.6,0.5,0.3],...
    'style','CheckBox',...
    'string','remove <',...
    'val',0);
HIO.labOpt_excluded_lessEdit = uicontrol('parent',HIO.labOpt_excluded,...
    'unit','norm',...
    'pos',[0.7,0.6,0.2,0.3],...
    'style','edit',...
    'enable','off');
HIO.labOpt_excluded_moreOpt = uicontrol('parent',HIO.labOpt_excluded,...
    'unit','norm',...
    'pos',[0.1,0.1,0.5,0.3],...
    'style','CheckBox',...
    'string','remove >',...
    'val',0);
HIO.labOpt_excluded_moreEdit = uicontrol('parent',HIO.labOpt_excluded,...
    'unit','norm',...
    'pos',[0.7,0.1,0.2,0.3],...
    'style','edit',...
    'enable','off');

HIO.labOpt_levelset = uibuttongroup('parent',HIO.labOpt,...
    'unit','norm',...
    'pos',[0.05,0.1,0.9,0.5],...
    'title','Level SetUp');
HIO.labOpt_level_discretevariable = uicontrol('parent',HIO.labOpt_levelset,...
    'unit','norm',...
    'pos',[0.1,0.9,0.8,0.08],...
    'style','rad',...
    'string','Discrete Variable');
HIO.labOpt_level_continuousvariable = uicontrol('parent',HIO.labOpt_levelset,...
    'unit','norm',...
    'pos',[0.1,0.8,0.8,0.08],...
    'style','rad',...
    'string','Continous Variable');
HIO.labOpt_level_ConVarBG = uibuttongroup('parent',HIO.labOpt_levelset,...
    'unit','norm',...
    'pos',[0.1,0.05,0.8,0.7],...
    'title','Continous Variable threshold');
HIO.labOpt_level_ConVar_middle = uicontrol('parent',HIO.labOpt_level_ConVarBG,...
    'unit','norm',...
    'pos',[0,2/3,0.5,1/3],...
    'style','rad',...
    'string','Mean');
HIO.labOpt_level_ConVar_median = uicontrol('parent',HIO.labOpt_level_ConVarBG,...
    'unit','norm',...
    'pos',[0.5,2/3,0.5,1/3],...
    'style','rad',...
    'string','Median');
HIO.labOpt_level_ConVar_quartile = uicontrol('parent',HIO.labOpt_level_ConVarBG,...
    'unit','norm',...
    'pos',[0,1/3,0.5,1/3],...
    'style','rad',...
    'string','quartile');
HIO.labOpt_level_ConVar_triPart = uicontrol('parent',HIO.labOpt_level_ConVarBG,...
    'unit','norm',...
    'pos',[0.5,1/3,0.5,1/3],...
    'style','rad',...
    'string','triPart');
HIO.labOpt_level_ConVar_freethreshold = uicontrol('parent',HIO.labOpt_level_ConVarBG,...
    'unit','norm',...
    'pos',[0,0/3,0.5,1/3],...
    'style','rad',...
    'string','Free Threshold');
HIO.labOpt_level_ConVar_freethresholdEdit = uicontrol('parent',HIO.labOpt_level_ConVarBG,...
    'unit','norm',...
    'pos',[0.5,0/3,0.5,1/3],...
    'style','edit');
HIO.labOpt_save = uicontrol('parent',HIO.labOpt,...
    'unit','norm',...
    'pos',[0.2,0,0.6,0.1],...
    'style','pushbutton',...
    'string','save Option');

%%
HIO.Showpbg = uibuttongroup('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.65,0.2,0.25,0.575],...
    'title','Label Showing');

HIO.lablistshow = uicontrol('parent',HIO.Showpbg,...
    'unit','norm',...
    'pos',[0.2,0.9,0.6,0.08],...
    'style','pushbutton',...
    'string','Show');
HIO.Showlab = axes('parent',HIO.Showpbg,...
    'unit','norm',...
    'pos',[0.05,0.05,0.9,0.8]);
axis(HIO.Showlab,'off')


%
HIO.spacesizet = uicontrol('parent',HIO.fig,...
     'unit','norm',...
     'pos',[0.2,0.05,0.1,0.05],...
     'style','text',...
     'string','Space Resolution');
HIO.spacesizee = uicontrol('parent',HIO.fig,...
     'unit','norm',...
     'pos',[0.3,0.05,0.1,0.05],...
     'style','edit',...
     'string','1');

HIO.OK = uicontrol('parent',HIO.fig,...
    'unit','norm',...
    'pos',[0.6,0.05,0.2,0.05],...
    'style','pushbutton',...
    'string','Finished!');

%%

set(HIO.Normsel,'val',1);

set(HIO.labOpt_level_ConVar_middle,'enable','off');
set(HIO.labOpt_level_ConVar_median,'enable','off');
set(HIO.labOpt_level_ConVar_quartile,'enable','off');
set(HIO.labOpt_level_ConVar_triPart,'enable','off');
set(HIO.labOpt_level_ConVar_freethreshold,'enable','off');
set(HIO.labOpt_level_ConVar_freethresholdEdit,'enable','off');

set(HIO.labOpt_del,'enable','off');
set(HIO.Normtype1,'enable','on')
set(HIO.Normtype2,'enable','on')
set(HIO.Normtype3,'enable','on')
set(HIO.Normtype4,'enable','on')
set(HIO.Normtype5,'enable','on')


set(HIO.labOpt_excluded_lessOpt,'callback',{@LowExcluded,HIO});
set(HIO.labOpt_excluded_moreOpt,'callback',{@HighExcluded,HIO});
set(HIO.labOpt_add,'callback',{@LabADD,HIO});
set(HIO.labOpt_del,'callback',{@LabDEL,HIO});
set(HIO.Normsel,'callback',{@HIOnormselfun,HIO});
set(HIO.inpb,'callback',{@HIOinput,HIO});
set(HIO.labOpt_level_discretevariable,'callback',{@selDisVar,HIO});
set(HIO.labOpt_level_continuousvariable,'callback',{@selConVar,HIO});
set(HIO.outpb,'callback',{@HIOoutput,HIO});
set(HIO.lablistshow,'callback',{@Showlist,HIO});
set(HIO.labOpt_save,'callback',{@SaveOpt,HIO});
set(HIO.lablist,'callback',{@loadopt,HIO});

set(HIO.OK,'callback',{@FinishJob,HIO,H});
end

function HIOnormselfun(varargin)
HIO = varargin{3};
val1 = get(HIO.Normsel,'val');
if val1==1
    set(HIO.Normtype1,'enable','on');
    set(HIO.Normtype2,'enable','on');
    set(HIO.Normtype3,'enable','on');
    set(HIO.Normtype4,'enable','on');
    set(HIO.Normtype5,'enable','on');
else    
    set(HIO.Normtype1,'enable','off');
    set(HIO.Normtype2,'enable','off');
    set(HIO.Normtype3,'enable','off');
    set(HIO.Normtype4,'enable','off');
    set(HIO.Normtype5,'enable','off');
end
end
function HIOinput(varargin)
HIO = varargin{3};
pgs = uigetdir(pwd,'InputDirectory');
set(HIO.inedit,'string',pgs);
HIO.para.Input = pgs;
end
function HIOoutput(varargin)
HIO = varargin{3};
pgs = uigetdir(pwd,'outputDirectory');
set(HIO.outedit,'string',pgs);
HIO.para.Output = pgs;
end
function LabADD(varargin)
HIO = varargin{3};
[answ] = inputdlg('Label name','Label name',1,{''});
[files pths] = uigetfile('*.txt',['Labels: ',answ{1}]);
StringofList = get(HIO.lablist,'string');
if isempty(StringofList)
    NumbofList = 1;
    StringofList{NumbofList,1} = ['Lab',sprintf('%03d',NumbofList),': ',answ{1}];   
else
    NumbofList = size(StringofList,1)+1;
    StringofList{NumbofList,1} = ['Lab',sprintf('%03d',NumbofList),': ',answ{1}];    
end
if ~isempty(dir([HIO.IOsetupdir,'LabList.mat']))
    load([HIO.IOsetupdir,'LabList.mat']);
end
LabList{NumbofList,1} = answ{1};
LabList{NumbofList,2} = fullfile(pths,files);
save([HIO.IOsetupdir,'LabList.mat'],'LabList');
set(HIO.lablist,'string',StringofList);
set(HIO.labOpt_del,'enable','on');

end
function LabDEL(varargin)
HIO = varargin{3};
% StringofList = get(HIO.lablist,'string');
Vals = get(HIO.lablist,'val');

load([HIO.IOsetupdir,'LabList.mat']);
nums = size(LabList,1);
% LabList = HIO.para.LabList;
%
LabList(Vals,:) =  [];
if isempty(LabList)
    set(HIO.labOpt_del,'enable','off');
    set(HIO.lablist,'string',[]);
    delete([HIO.IOsetupdir,'LabList.mat']);
    delete([HIO.IOsetupdir,'Var',num2str(Vals),'Opt.mat'])
else
    for i = 1:size(LabList,1)
        StringofList{i,1}= ['Lab',sprintf('%03d',i),': ',LabList{i,1}];
    end
    
    save([HIO.IOsetupdir,'LabList.mat'],'LabList')
    set(HIO.lablist,'val',1);
    set(HIO.lablist,'string',StringofList);
    delete([HIO.IOsetupdir,'Var',num2str(Vals),'Opt.mat']);
    for i = Vals:nums-1
        copyfile([HIO.IOsetupdir,'Var',num2str(i+1),'Opt.mat'],[HIO.IOsetupdir,'Var',num2str(i),'Opt.mat']);
        delete([HIO.IOsetupdir,'Var',num2str(i+1),'Opt.mat']);        
    end
end
end
function LowExcluded(varargin)
HIO = varargin{3};
VAL = get(HIO.labOpt_excluded_lessOpt,'val');
if VAL==1
    set(HIO.labOpt_excluded_lessEdit,'enable','on');
else
    set(HIO.labOpt_excluded_lessEdit,'enable','off');
end
end
function HighExcluded(varargin)
HIO = varargin{3};
VAL = get(HIO.labOpt_excluded_moreOpt,'val');
if VAL==1
    set(HIO.labOpt_excluded_moreEdit,'enable','on');
else
    set(HIO.labOpt_excluded_moreEdit,'enable','off');
end
end
function selDisVar(varargin)
HIO = varargin{3};
vals = get(HIO.labOpt_level_discretevariable,'val');
if vals==1
    set(HIO.labOpt_level_ConVar_middle,'enable','off');
    set(HIO.labOpt_level_ConVar_median,'enable','off');
    set(HIO.labOpt_level_ConVar_quartile,'enable','off');
    set(HIO.labOpt_level_ConVar_triPart,'enable','off');
    set(HIO.labOpt_level_ConVar_freethreshold,'enable','off');
    set(HIO.labOpt_level_ConVar_freethresholdEdit,'enable','off');
end
end
function selConVar(varargin)
HIO = varargin{3};
vals = get(HIO.labOpt_level_continuousvariable,'val');

if vals==1
    set(HIO.labOpt_level_ConVar_middle,'enable','on');
    set(HIO.labOpt_level_ConVar_median,'enable','on');
    set(HIO.labOpt_level_ConVar_quartile,'enable','on');
    set(HIO.labOpt_level_ConVar_triPart,'enable','on');
    set(HIO.labOpt_level_ConVar_freethreshold,'enable','on');
    set(HIO.labOpt_level_ConVar_freethresholdEdit,'enable','on');
end
end
function Showlist(varargin)
HIO = varargin{3};
load([HIO.IOsetupdir,'LabList.mat'])
val = get(HIO.lablist,'val');
% StringofList = get(HIO.lablist,'string');
% LabList = HIO.para.LabList;
ExcLab1 = get(HIO.labOpt_excluded_lessOpt,'val');
ExcLab2 = get(HIO.labOpt_excluded_moreOpt,'val');
ExcLab1_valtemp = get(HIO.labOpt_excluded_lessEdit,'string');
ExcLab2_valtemp = get(HIO.labOpt_excluded_moreEdit,'string');
if ExcLab1==1
    para.Exclude{val,1} = str2num(ExcLab1_valtemp);
else
    para.Exclude{val,1} = -inf;
end
if ExcLab2==1
    para.Exclude{val,2} = str2num(ExcLab2_valtemp);
else
    para.Exclude{val,2} = inf;
end
Dat = load(LabList{val,2});

% Dat2 = Dat;

ind1 = find(Dat<para.Exclude{val,2}&Dat>para.Exclude{val,1});
ind2 = 1:length(Dat);
ind2(ind1) = [];

DISVARlab = get(HIO.labOpt_level_discretevariable,'val');
CONVARlab = get(HIO.labOpt_level_continuousvariable,'val');
CONVARtype1 = get(HIO.labOpt_level_ConVar_middle,'val');
CONVARtype2 = get(HIO.labOpt_level_ConVar_median,'val');
CONVARtype3 = get(HIO.labOpt_level_ConVar_quartile,'val');
CONVARtype4 = get(HIO.labOpt_level_ConVar_triPart,'val');
CONVARtype5 = get(HIO.labOpt_level_ConVar_freethreshold,'val');
CONVARtype5dattemp = get(HIO.labOpt_level_ConVar_freethresholdEdit,'string');
CONVARtype5val = str2num(CONVARtype5dattemp);

if DISVARlab
    Dat2 = Dat(ind1);
    maxdat = max(Dat2);
    Dat2 = (Dat2-min(Dat2))/maxdat*0.9+0.05;
    [Dat2re Dat2ind] = sort(Dat2);
    DATshow(ind1,1) = Dat2;
    DATshow(ind2,1) = 2;
    DATshow(1:length(ind1),2) = Dat2re;
    if length(ind1)<length(Dat)
        DATshow(length(ind1)+1:length(Dat),2) = 2;
    end
    image(HIO.Showlab,DATshow,'CDataMapping','scaled');
    colormap(HIO.Showlab,HIO.colorforlab)
    set(HIO.Showlab,'CLim',[0,2])
    set(HIO.Showpbg,'title',['Dis Var: ',LabList{val,1}])
    
else
    Dat2 = Dat(ind1);
    maxdat = max(Dat2);
    Dat2 = (Dat2-min(Dat2))/(maxdat-min(Dat2))*0.9+0.05;
    [Dat2re Dat2ind] = sort(Dat2);
    DATshow(ind1,1) = Dat2;
    DATshow(ind2,1) = 2;
    DATshow(1:length(ind1),2) = Dat2re;
    if length(ind1)<length(Dat)
        DATshow(length(ind1)+1:length(Dat),2) = 2;
    end
    
    DATA = Dat(ind1);
    if CONVARtype1
        midv = mean(DATA);
        thr = midv;
    end
    if CONVARtype2
        medv = median(DATA);
        thr = medv;
    end
    if CONVARtype3
        totalnum = length(ind1);
        medv = median(DATA);
        [dat1re dat1ord] = sort(DATA);
        firstsq = round(totalnum/4);
        thirdsq = totalnum-firstsq;
        thr(1) = dat1re(firstsq);
        thr(2) = medv;
        thr(3) = dat1re(thirdsq);
    end
    if CONVARtype4        
        totalnum = length(ind1);
        [dat1re dat1ord] = sort(DATA);
        firsttri = round(totalnum/3);
        secondtri = round(totalnum/3*2);
        thr(1) = dat1re(firsttri);
        thr(2) = dat1re(secondtri);
    end
    if CONVARtype5
        thr = CONVARtype5val;
    end
    thrnum = length(thr);
    thrU(1) = min(DATA)-1;
    thrU(2:thrnum+1) = thr;
    thrU(thrnum+2) = max(DATA);
    for i = 1:thrnum+1
        INDTHR{i} = find(DATA<=thrU(i+1)&DATA>thrU(i));
        OutshowT(ind1(INDTHR{i})) = (i-0.5)/(thrnum+1);
    end
    if ~isempty(ind2)
        OutshowT(ind2) = 2;
    end
    [OutshowTs,OutshowTsind] = sort(OutshowT);
    DATshow(:,3) = OutshowT;
    DATshow(:,4) = OutshowTs;
    
    image(HIO.Showlab,DATshow,'CDataMapping','scaled');
    colormap(HIO.Showlab,HIO.colorforlab)
    set(HIO.Showlab,'CLim',[0,2])
    set(HIO.Showpbg,'title',['Dis Var: ',LabList{val,1}])
end


end
function SaveOpt(varargin)
HIO = varargin{3};
%%
load([HIO.IOsetupdir,'LabList.mat'])
val = get(HIO.lablist,'val');
% StringofList = get(HIO.lablist,'string');
% LabList = HIO.para.LabList;
ExcLab1 = get(HIO.labOpt_excluded_lessOpt,'val');
ExcLab2 = get(HIO.labOpt_excluded_moreOpt,'val');
ExcLab1_valtemp = get(HIO.labOpt_excluded_lessEdit,'string');
ExcLab2_valtemp = get(HIO.labOpt_excluded_moreEdit,'string');
if ExcLab1==1
    para.Exclude(1,1) = str2num(ExcLab1_valtemp);
else
    para.Exclude(1,1) = -inf;
    set(HIO.labOpt_excluded_lessEdit,'string',[]);
end
if ExcLab2==1
    para.Exclude(1,2) = str2num(ExcLab2_valtemp);
else
    para.Exclude(1,2) = inf;
    set(HIO.labOpt_excluded_moreEdit,'string',[]);
end
Dat = load(LabList{val,2});

% Dat2 = Dat;

ind1 = find(Dat<para.Exclude(1,2)&Dat>para.Exclude(1,1));
ind2 = 1:length(Dat);
ind2(ind1) = [];

DISVARlab = get(HIO.labOpt_level_discretevariable,'val');
CONVARlab = get(HIO.labOpt_level_continuousvariable,'val');
CONVARtype1 = get(HIO.labOpt_level_ConVar_middle,'val');
CONVARtype2 = get(HIO.labOpt_level_ConVar_median,'val');
CONVARtype3 = get(HIO.labOpt_level_ConVar_quartile,'val');
CONVARtype4 = get(HIO.labOpt_level_ConVar_triPart,'val');
CONVARtype5 = get(HIO.labOpt_level_ConVar_freethreshold,'val');
CONVARtype5dattemp = get(HIO.labOpt_level_ConVar_freethresholdEdit,'string');
CONVARtype5val = str2num(CONVARtype5dattemp);

if DISVARlab
    DATA = Dat(ind1);
    thr = unique(DATA);
    Optsave.dattype = DISVARlab;
    Optsave.thr = thr;
    for i = 1:length(thr)
        datalab{i} = ind1(DATA==thr(i));
    end
    Optsave.datalab = datalab;
    if ~isempty(ind2)
        Optsave.Excludedlab = ind2;
    end
    Optsave.para = para;
else
    
    Optsave.dattype = DISVARlab;
    
    Dat2 = Dat(ind1);
    maxdat = max(Dat2);
    Dat2 = (Dat2-min(Dat2))/maxdat*0.9+0.05;
    [Dat2re Dat2ind] = sort(Dat2);
    DATshow(ind1,1) = Dat2;
    DATshow(ind2,1) = 2;
    DATshow(1:length(ind1),2) = Dat2re;
    if length(ind1)<length(Dat)
        DATshow(length(ind1)+1:length(Dat),2) = 2;
    end
    
    DATA = Dat(ind1);
    if CONVARtype1
        midv = mean(DATA);
        thr = midv;
        Optsave.con_type = 1;
    end
    if CONVARtype2
        medv = median(DATA);
        thr = medv;
        
        Optsave.con_type = 2;
    end
    if CONVARtype3
        totalnum = length(ind1);
        medv = median(DATA);
        [dat1re dat1ord] = sort(DATA);
        firstsq = round(totalnum/4);
        thirdsq = totalnum-firstsq;
        thr(1) = dat1re(firstsq);
        thr(2) = medv;
        thr(3) = dat1re(thirdsq);
        
        Optsave.con_type = 3;
    end
    if CONVARtype4        
        totalnum = length(ind1);
        [dat1re dat1ord] = sort(DATA);
        firsttri = round(totalnum/3);
        secondtri = round(totalnum/3*2);
        thr(1) = dat1re(firsttri);
        thr(2) = dat1re(secondtri);
        
        Optsave.con_type = 4;
    end
    if CONVARtype5
        thr = CONVARtype5val;
        
        Optsave.con_type = 5;
        Optsave.con_input = thr;
    end
    Optsave.thr = thr;
    thrnum = length(thr);
    thrU(1) = min(DATA)-1;
    thrU(2:thrnum+1) = thr;
    thrU(thrnum+2) = max(DATA);
    for i = 1:thrnum+1
        datalab{i} = ind1(find(DATA<=thrU(i+1)&DATA>thrU(i)));        
    end
    if ~isempty(ind2)
        Optsave.Excludedlab = ind2;
    end
    Optsave.datalab = datalab;
    Optsave.para = para;
end
save([HIO.IOsetupdir,'Var',num2str(val),'Opt.mat'],'Optsave');
end
function loadopt(varargin)
HIO = varargin{3};

load([HIO.IOsetupdir,'LabList.mat'])
val = get(HIO.lablist,'val');
if isempty(dir([HIO.IOsetupdir,'Var',num2str(val),'Opt.mat']))
    uiwait(msgbox('Please confirm the options, and click [save Option]'))
    
else
    load([HIO.IOsetupdir,'Var',num2str(val),'Opt.mat'])
    para = Optsave.para;
    if para.Exclude(1,1)>-inf
        set(HIO.labOpt_excluded_lessOpt,'val',1);
        set(HIO.labOpt_excluded_lessEdit,'string',num2str(para.Exclude(1,1)));
    else
        set(HIO.labOpt_excluded_lessOpt,'val',0);
    end
    if para.Exclude(1,2)<inf
        set(HIO.labOpt_excluded_moreOpt,'val',1);
        set(HIO.labOpt_excluded_moreEdit,'string',num2str(para.Exclude(1,2)));
    else        
        set(HIO.labOpt_excluded_moreOpt,'val',0);
    end
    if Optsave.dattype % 1 for dis , 0  for cont
        set(HIO.labOpt_level_continuousvariable,'val',0);
        set(HIO.labOpt_level_ConVar_freethreshold,'enable','off');
        set(HIO.labOpt_level_ConVar_freethresholdEdit,'enable','off');
        set(HIO.labOpt_level_ConVar_median,'enable','off');
        set(HIO.labOpt_level_ConVar_middle,'enable','off');
        set(HIO.labOpt_level_ConVar_quartile,'enable','off');
        set(HIO.labOpt_level_ConVar_triPart,'enable','off');
        set(HIO.labOpt_level_discretevariable,'val',1);
    else        
        set(HIO.labOpt_level_continuousvariable,'val',1);
        set(HIO.labOpt_level_ConVar_freethreshold,'enable','on');
        set(HIO.labOpt_level_ConVar_freethresholdEdit,'enable','on');
        set(HIO.labOpt_level_ConVar_median,'enable','on');
        set(HIO.labOpt_level_ConVar_middle,'enable','on');
        set(HIO.labOpt_level_ConVar_quartile,'enable','on');
        set(HIO.labOpt_level_ConVar_triPart,'enable','on');
        set(HIO.labOpt_level_discretevariable,'val',0);
        switch Optsave.con_type
            case 1
                set(HIO.labOpt_level_ConVar_middle,'val',1);
                set(HIO.labOpt_level_ConVar_median,'val',0);
                set(HIO.labOpt_level_ConVar_quartile,'val',0);
                set(HIO.labOpt_level_ConVar_triPart,'val',0);
                set(HIO.labOpt_level_ConVar_freethreshold,'val',0);
                set(HIO.labOpt_level_ConVar_freethresholdEdit,'string','');
            case 2
                set(HIO.labOpt_level_ConVar_middle,'val',0);
                set(HIO.labOpt_level_ConVar_median,'val',1);
                set(HIO.labOpt_level_ConVar_quartile,'val',0);
                set(HIO.labOpt_level_ConVar_triPart,'val',0);
                set(HIO.labOpt_level_ConVar_freethreshold,'val',0);
                set(HIO.labOpt_level_ConVar_freethresholdEdit,'string','');
                
            case 3
                set(HIO.labOpt_level_ConVar_middle,'val',0);
                set(HIO.labOpt_level_ConVar_median,'val',0);
                set(HIO.labOpt_level_ConVar_quartile,'val',1);
                set(HIO.labOpt_level_ConVar_triPart,'val',0);
                set(HIO.labOpt_level_ConVar_freethreshold,'val',0);
                set(HIO.labOpt_level_ConVar_freethresholdEdit,'string','');
                
            case 4
                set(HIO.labOpt_level_ConVar_middle,'val',0);
                set(HIO.labOpt_level_ConVar_median,'val',0);
                set(HIO.labOpt_level_ConVar_quartile,'val',0);
                set(HIO.labOpt_level_ConVar_triPart,'val',1);
                set(HIO.labOpt_level_ConVar_freethreshold,'val',0);
                set(HIO.labOpt_level_ConVar_freethresholdEdit,'string','');
                
            case 5
                set(HIO.labOpt_level_ConVar_middle,'val',0);
                set(HIO.labOpt_level_ConVar_median,'val',0);
                set(HIO.labOpt_level_ConVar_quartile,'val',0);
                set(HIO.labOpt_level_ConVar_triPart,'val',0);
                set(HIO.labOpt_level_ConVar_freethreshold,'val',1);
                set(HIO.labOpt_level_ConVar_freethresholdEdit,'string',mat2str(Optsave.con_input));
        end
        
    end
    
end
end
function FinishJob(varargin)
HIO = varargin{3};
H = varargin{4};
load([HIO.IOsetupdir,'LabList.mat']);
NoVar = size(LabList,1);
for i = 1:NoVar
    load([HIO.IOsetupdir,'Var',num2str(i),'Opt.mat'])
    ParaOut(i).dattype = Optsave.dattype;
    ParaOut(i).datalab = Optsave.datalab;
    if isfield(Optsave,'Excludedlab')
        ParaOut(i).Excludedlab = Optsave.Excludedlab;
    else
        ParaOut(i).Excludedlab = [];        
    end
    ParaOut(i).thr = Optsave.thr;
    ParaOut(i).LabName = LabList{i,1};
    ParaOut(i).Labtxt = LabList{i,2};
end
SetUpPara.Indir = get(HIO.inedit,'string');
SetUpPara.Outdir = get(HIO.outedit,'string');
SetUpPara.ParaOut = ParaOut;
Vsize = get(HIO.spacesizee,'string');
SetUpPara.OutVsize = str2num(Vsize);
SetUpPara.NormOpt = get(HIO.Normsel,'val');
SetUpPara.Normtype1 = get(HIO.Normtype1,'val');
SetUpPara.Normtype2 = get(HIO.Normtype2,'val');
SetUpPara.Normtype3 = get(HIO.Normtype3,'val');
SetUpPara.Normtype4 = get(HIO.Normtype4,'val');
SetUpPara.Normtype5 = get(HIO.Normtype5,'val');
save([SetUpPara.Outdir,filesep,'SetUpPara.mat'],'SetUpPara');
save([HIO.IOsetupdir,filesep,'SetUpPara.mat'],'SetUpPara');
cd(SetUpPara.Outdir);
close(HIO.fig);
movegui(H.fig,'center');
end