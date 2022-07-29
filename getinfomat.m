function Pvalp = getinfomat(cutnum,permclust_dir,inds,iperm)

Pvalp = zeros(length(inds),1);
for icut = 1:cutnum
    load([permclust_dir,'Piece',sprintf('%05d',icut),'.mat']);
    INDS = CUTpiece{icut};
    Pvalp(INDS) = Pvalperm(:,iperm);
end
end