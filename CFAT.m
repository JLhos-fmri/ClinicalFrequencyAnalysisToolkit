function CFAT
% clinic frequency analysis main function
% Version 1: fisrt build in 02 Dec. 2021
% Updata on 11 Mar. 2022
% Qiang Xu, 
% Department of Radiology, Jinling hospital
% Medical School of Nanjing University
% fmrixuq@126.com

%%=========================================================================
% part 1: I/O, input & setup labels, standard space setup
% part 2: template setup
% part 3: Check and Run (statistical analysis included)
% part 4: Show Results
% part 5: Utilities (regroup dataset)
% part 6: Exit
%%=========================================================================
clc;
clear all;
close all;
uiwait(msgbox({'Clinical Frequency Analysis Toolkit',...
    'Developed by Department of Radiology, Jinling hospital, Medical School, Nanjing University',...
    'If you have any comments, please contract fmrixuq@126.com'}));
disp('Clinical Frequency Analysis Toolkit!')
disp('Developed by Department of Radiology, Jinling hospital, Medical School, Nanjing University')
disp('If you have any comments, please contract fmrixuq@126.com')
H.fig = figure('unit','norm',...
    'pos',[0.3,0.3,0.4,0.4],...
    'menubar','none',...
    'color','w',...
    'name','Clinical Frequency Analysis Toolkit in Nervous System');

H.part1pb = uicontrol('parent',H.fig,...
    'unit','norm',...
    'pos',[0.1,0.65,0.2,0.2],...
    'style','pushbutton',...
    'string','I/O & Setup');

H.part2pb = uicontrol('parent',H.fig,...
    'unit','norm',...
    'pos',[0.4,0.65,0.2,0.2],...
    'style','pushbutton',...
    'string','Template Select',...
    'enable','off');

H.part3pb = uicontrol('parent',H.fig,...
    'unit','norm',...
    'pos',[0.7,0.65,0.2,0.2],...
    'style','pushbutton',...
    'string','Check & Run',...
    'enable','off');

H.part4pb = uicontrol('parent',H.fig,...
    'unit','norm',...
    'pos',[0.1,0.15,0.2,0.2],...
    'style','pushbutton',...
    'string','Make Result for Show');

H.part5pb = uicontrol('parent',H.fig,...
    'unit','norm',...
    'pos',[0.4,0.15,0.2,0.2],...
    'style','popupmenu',...
    'string',{'Utilties','Regroup','Check ROIs','batch mode'});

H.part6pb = uicontrol('parent',H.fig,...
    'unit','norm',...
    'pos',[0.7,0.15,0.2,0.2],...
    'style','pushbutton',...
    'string','Exit');

set(H.part1pb,'callback',{@IOselAndSetUp,H});
set(H.part2pb,'callback',{@SelectTemplate,H});
set(H.part3pb,'callback',{@CheckAndRun,H});
set(H.part4pb,'callback',{@Showmode,H});
set(H.part5pb,'callback',{@UtiltiesMode,H});
set(H.part6pb,'callback',{@ExitFun,H});
end
function IOselAndSetUp(varargin)
H = varargin{3};
set(H.part2pb,'enable','on');
CFA_inout(H);
end
function SelectTemplate(varargin)
H = varargin{3};
set(H.part3pb,'enable','on');
CFA_template(H);
end
function CheckAndRun(varargin)
H = varargin{3};
CFA_Checkandrun(H);
end
function Showmode(varargin)
H = varargin{3};
CFA_MakeRes;
end
function UtiltiesMode(varargin)
H = varargin{3};
vals = get(H.part5pb,'val');
switch vals
    case 2 % ReGroup
        CFA_ReGroupData
    case 3 % CheckROIs
        CFA_CheckData
    case 4 % batch
        BatchModeShow
    otherwise
        
end
end
function ExitFun(varargin)
H = varargin{3};
close(H.fig);

end