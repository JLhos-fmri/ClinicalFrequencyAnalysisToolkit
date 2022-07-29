function CFA_SingleGroupPval_Perm(DATtemp,VoxelDir,LabName,V,j,IDtemp)
subnum = size(DATtemp,2);
Indexist = find(sum(DATtemp,2));
% PermDat = DATtemp(:,Indexist);
% full(sum(DATtemp,2))
bs_num = 5000;
bs_mat = ceil(rand(bs_num,subnum)*subnum);
% bs_mat(bs_mat==0) = 1;
% Edat = zeros(length(Indexist),bs_num);

Seplength = 50000;
sepnum = floor(length(Indexist)/Seplength);
if sepnum>0
    for i = 1:sepnum
        Sepind{i,1} = 1+Seplength*(i-1):Seplength*i;
    end
end
if Indexist>Seplength*sepnum
    sepnum = sepnum+1;
    Sepind{sepnum,1} = Seplength*(sepnum-1)+1:Indexist;
end
  
Pvals = zeros(length(Indexist),1);
Pvals2 = zeros(length(Indexist),1);
tic
for isep = 1:sepnum
    isep
    voxnum = Sepind{isep,1};
    parfor i = 1:bs_num
        Edat(:,i) = full(sum(DATtemp(Indexist(voxnum),bs_mat(i,:)),2));
    end
    Origdat = full(sum(DATtemp(Indexist(voxnum),:),2));
    parfor i = 1:length(voxnum)
        [mu sig muci sigci] = normfit(Edat(i,:));
        Pval(i,1) = normcdf(Origdat(i,1),mu,sig);
    end
    parfor i = 1:length(voxnum)
        [Freqv Centv] = hist(Edat(i,:));
        freqpercent = cumsum(Freqv)/sum(Freqv);
        indss = min(find(Origdat(i,1)>Centv));
        Pval2(i,1) = freqpercent(indss);
    end
    Pvals(voxnum) = Pval;
    Pvals2(voxnum) = Pval2;
    clear Pval Edat Pval2
    toc
end
OutPval = zeros(V(1).dim);
OutPval(Indexist) = Pvals;
DynamicBC_write_NIFTI(OutPval,V(1),[VoxelDir,filesep,'Pval_',LabName,'_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.nii'])
OutPval = zeros(V(1).dim);
OutPval(Indexist) = Pvals2;
DynamicBC_write_NIFTI(OutPval,V(1),[VoxelDir,filesep,'Pval_bs_',LabName,'_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.nii'])
% [VoxelDir,filesep,LabName,'_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.nii']
end