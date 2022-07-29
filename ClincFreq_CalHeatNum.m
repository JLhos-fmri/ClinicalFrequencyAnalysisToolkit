function ClincFreq_CalHeatNum(INDIR,OutHeatdir,SetUpPara,Para,AtlasTempU,dbmask)
[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(INDIR);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
D = D>0;

for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    HeatDir{i} = [OutHeatdir,filesep,LabName,filesep];
    mkdir(HeatDir{i});
    if Para.Heat.V
        VoxelDir = [HeatDir{i},'VoxelWise',filesep];
        mkdir(VoxelDir);
    end
    if Para.Heat.R
        ROIDir = [HeatDir{i},'ROIWise',filesep];
        mkdir(ROIDir);
        for iROI = 1:size(AtlasTempU,1)
            mkdir([ROIDir,AtlasTempU{iROI,3},filesep]);
        end
    end
    for j = 1:length(SetUpPara.ParaOut(i).datalab)
        IDtemp = SetUpPara.ParaOut(i).datalab{j};
        GroupNum(j,1) = length(IDtemp);
        DATtemp = D(:,IDtemp);
        
        if Para.Heat.V
            DATtempV = reshape(full(sum(DATtemp,2)),V(1).dim).*dbmask;
            DynamicBC_write_NIFTI(DATtempV,V(1),[VoxelDir,filesep,'Freq_',LabName,'_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.nii'])
%             CFA_SingleGroupPval_Perm(DATtemp,VoxelDir,LabName,V(1),j,IDtemp);
        end
        if Para.Heat.R
            for iROI = 1:size(AtlasTempU,1)
                PercentVol = Para.Heat.Rper/100;
                indtab = AtlasTempU{iROI,6};
                for k = 1:length(IDtemp)
                    DATsingle = full(DATtemp(:,k));
                    for l = 1:length(indtab)
                        % Marked by the size of Lesion mask.  Affected
                        % Voxel number of ROI divided the whole number of
                        % Lesion.
                        INDinvolv1(k,l) = nnz(DATsingle(indtab{l}))/nnz(DATsingle);
                        % Marked by the size of Atlas. Affected Voxel
                        % number of ROI divided the voxel number of Atlas.
                        INDinvolv2(k,l) = nnz(DATsingle(indtab{l}))/length(indtab{l});
                        % 
                        INDinvolv3(k,l) = max(INDinvolv1(k,l),INDinvolv2(k,l));
                    end
                end
                IND_involve1 = INDinvolv1>PercentVol;
                IND_involve2 = INDinvolv2>PercentVol;
                IND_involve3 = INDinvolv3>PercentVol;
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedROI.csv'],IND_involve1);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedAtlas.csv'],IND_involve2);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedBoth.csv'],IND_involve3);
                NamelistSave = namelist(IDtemp);
%                 DATe = zeros(V(1).dim);
                save([ROIDir,AtlasTempU{iROI,3},filesep,'Single_Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'.mat'],'NamelistSave','V','IND_involve1','IND_involve2','IND_involve3');
                clear NamtlistSave
                
                
                IND_total1 = sum(IND_involve1,1);
                IND_total2 = sum(IND_involve2,1);
                IND_total3 = sum(IND_involve3,1);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedROI.csv'],IND_total1);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedAtlas.csv'],IND_total2);
                csvwrite([ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-total',num2str(length(IDtemp)),'-markedBoth.csv'],IND_total3);
                
                
                DATu = zeros(V(1).dim);
                for l = 1:length(indtab)
                    DATu(indtab{l}) = IND_total1(l);
                end
                DynamicBC_write_NIFTI(DATu,V(1),[ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-ROImarked.nii']);
                DATu = zeros(V(1).dim);
                for l = 1:length(indtab)
                    DATu(indtab{l}) = IND_total2(l);
                end
                DynamicBC_write_NIFTI(DATu,V(1),[ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-Atlasmarked.nii']);
                DATu = zeros(V(1).dim);
                for l = 1:length(indtab)
                    DATu(indtab{l}) = IND_total3(l);
                end
                DynamicBC_write_NIFTI(DATu,V(1),[ROIDir,AtlasTempU{iROI,3},filesep,'Group',sprintf('%03d',j),'-Bothmarked.nii']);
                
                clear INDinvolv1 INDinvolv2 INDinvolv3 IND_involve1 IND_involve2 IND_involve3 IND_total1 IND_total2 IND_total3 indtab
            end
            
        end
        
    end
    
end

end