function [PvalpermV0001C,PvalpermV0005C,PvalpermV001C,PvalpermV005C] = getclusterinfo(PvalpermV005,PvalpermV001,PvalpermV0005,PvalpermV0001,iperm);


Coninf = bwconncomp(PvalpermV005);
nums = Coninf.NumObjects;
if nums>0
    for i = 1:nums
        PvalpermV005C(i,1) = length(Coninf.PixelIdxList{i});
    end
else
    PvalpermV005C = 0;
end

Coninf = bwconncomp(PvalpermV001);
nums = Coninf.NumObjects;
if nums>0
    for i = 1:nums
        PvalpermV001C(i,1) = length(Coninf.PixelIdxList{i});
    end
else
    PvalpermV001C = 0;
end

Coninf = bwconncomp(PvalpermV0005);
nums = Coninf.NumObjects;
if nums>0
    for i = 1:nums
        PvalpermV0005C(i,1) = length(Coninf.PixelIdxList{i});
    end
else
    PvalpermV0005C = 0;
end
Coninf = bwconncomp(PvalpermV0001);
nums = Coninf.NumObjects;
if nums>0
    for i = 1:nums
        PvalpermV0001C(i,1) = length(Coninf.PixelIdxList{i});
    end
else
    PvalpermV0001C = 0;
end
%
% Coninf = bwconncomp(PvalpermV005);
% nums = Coninf.NumObjects;
% for i = 1:nums
%     PvalpermV005C(iperm,i) = length(Coninf.PixelIdxList{i});
% end
% Coninf = bwconncomp(PvalpermV001);
% nums = Coninf.NumObjects;
% for i = 1:nums
%     PvalpermV001C(iperm,i) = length(Coninf.PixelIdxList{i});
% end
% Coninf = bwconncomp(PvalpermV0005);
% nums = Coninf.NumObjects;
% for i = 1:nums
%     PvalpermV0005C(iperm,i) = length(Coninf.PixelIdxList{i});
% end
% Coninf = bwconncomp(PvalpermV0001);
% nums = Coninf.NumObjects;
% for i = 1:nums
%     PvalpermV0001C(iperm,i) = length(Coninf.PixelIdxList{i});
% end
end