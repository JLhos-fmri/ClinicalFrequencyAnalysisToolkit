function ClusterNum = PermClusterInfo(PvalpermV0001Ct)
for i = 1:length(PvalpermV0001Ct)
    Pvaltemps(i,1:length(PvalpermV0001Ct{i})) = PvalpermV0001Ct{i};
end
AAA = Pvaltemps(Pvaltemps>0);
BBB = unique(AAA);
for i = 1:length(BBB)
    noaaa(i) = length(find(AAA==BBB(i)));
end
CCC = cumsum(noaaa)/sum(noaaa);


inds = min(find(CCC>0.95));
if ~isempty(inds)
    ClusterNum.P005 = BBB(inds);
else
    ClusterNum.P005 = 0;
end

inds = min(find(CCC>0.99));
if ~isempty(inds)
    ClusterNum.P001 = BBB(inds);
else
    ClusterNum.P001 = 0;
end

inds = min(find(CCC>0.995));
if ~isempty(inds)
    ClusterNum.P0005 = BBB(inds);
else
    ClusterNum.P0005 = 0;
end

inds = min(find(CCC>0.999));
if ~isempty(inds)
    ClusterNum.P0001 = BBB(inds);
else
    ClusterNum.P0001 = 0;
end
end