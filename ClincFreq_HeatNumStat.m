function ClincFreq_HeatNumStat(INDIR,OutHeatdir,SetUpPara,Para,AtlasTempU)
[V,D,namelist] = Dynamic_read_dir_NIFTI_sparse(INDIR);
D(isnan(D)) = 0;
D(isinf(D)) = 0;
D = D>0;
for i = 1:length(SetUpPara.ParaOut)
    LabName = SetUpPara.ParaOut(i).LabName;
    HeatDir{i} = [OutHeatdir,filesep,LabName,filesep];
    % if only one group
    if Para.Heat.V&&length(dir([HeatDir{i},'VoxelWise',filesep,'Freq_*.nii']))<2||...
            Para.Heat.R&&length(dir([HeatDir{i},'ROIWise',filesep,AtlasTempU{1,3},filesep,'Group*-markedAtlas.csv']))<2
        continue
    end
    
    if Para.Heat.V
        VoxelDir = [HeatDir{i},'VoxelWise',filesep];
        %         mkdir(VoxelDir);
        if Para.Heat.StChisquare
            files = dir([VoxelDir,'Freq_*.nii']);
            
            VoxelDirChi2test = [HeatDir{i},'VoxelWise_Chi2test',filesep];
            mkdir(VoxelDirChi2test);
            
            for ifile = 1:length(files)
                filetemp = files(ifile).name;
                ind = find(filetemp=='-');
                GroupNum(ifile,1) = str2num(filetemp(ind(end)+6:end-4));
                [Vraw Draw] = Dynamic_read_dir_NIFTI([VoxelDir,files(ifile).name]);
                Dat(ifile,:) = Draw;
            end
            C = nchoosek(1:length(files),2);
            
            DATAs = sum(Dat);
            INDDATA = find(DATAs);
            DATAu = Dat(:,INDDATA);
            DATAu2 = repmat(GroupNum,1,length(INDDATA))-DATAu;
            p_ANOVA = zeros(length(INDDATA),1);
            Q_ANOVA = zeros(length(INDDATA),1);
            p_two = zeros(length(INDDATA),size(C,1));
            q_two = zeros(length(INDDATA),size(C,1));
            or_two = zeros(length(INDDATA),size(C,1));
            or_twolow = zeros(length(INDDATA),size(C,1));
            or_twohigh = zeros(length(INDDATA),size(C,1));
            parfor idat = 1:length(INDDATA)
                x = [DATAu(:,idat),DATAu2(:,idat)];
                [p_ANOVA(idat),Q_ANOVA(idat),p_two(idat,:),q_two(idat,:),or_two(idat,:),or_twolow(idat,:),or_twohigh(idat,:)] = ClincFreq_Chi2test(x,C)
            end
            p_ANOVA(isnan(p_ANOVA)) = 0;
            Q_ANOVA(isnan(Q_ANOVA)) = 0;
            p_two(isnan(p_two)) = 0;
            q_two(isnan(q_two)) = 0;
            or_two(isnan(or_two)) = 0;
            or_twolow(isnan(or_twolow)) = 0;
            or_twohigh(isnan(or_twohigh)) = 0;
            if size(C,1)>1
                DAT_panova = zeros(Vraw.dim);
                DAT_panova(INDDATA) = p_ANOVA;
                DAT_qanova = zeros(Vraw.dim);
                DAT_qanova(INDDATA) = Q_ANOVA;
                DynamicBC_write_NIFTI(DAT_panova,Vraw,[VoxelDirChi2test,'ANOVA_p.nii']);
                DynamicBC_write_NIFTI(DAT_qanova,Vraw,[VoxelDirChi2test,'ANOVA_q.nii']);
            end
            for iC = 1:size(C,1)
                dvaltemp = zeros(Vraw.dim);
                dvaltemp(INDDATA) = q_two(:,iC);
                DynamicBC_write_NIFTI(dvaltemp,Vraw,[VoxelDirChi2test,['Group',sprintf('%03d',C(iC,1)),'vsGroup',sprintf('%03d',C(iC,2))],'_q.nii']);
                dvaltemp = zeros(Vraw.dim);
                dvaltemp(INDDATA) = p_two(:,iC);
                DynamicBC_write_NIFTI(dvaltemp,Vraw,[VoxelDirChi2test,['Group',sprintf('%03d',C(iC,1)),'vsGroup',sprintf('%03d',C(iC,2))],'_p.nii']);
                
                dvaltemp = zeros(Vraw.dim);
                dvaltemp(INDDATA) = or_two(:,iC);
                DynamicBC_write_NIFTI(dvaltemp,Vraw,[VoxelDirChi2test,['Group',sprintf('%03d',C(iC,1)),'vsGroup',sprintf('%03d',C(iC,2))],'_OR.nii']);
                dvaltemp = zeros(Vraw.dim);
                dvaltemp(INDDATA) = or_twolow(:,iC);
                DynamicBC_write_NIFTI(dvaltemp,Vraw,[VoxelDirChi2test,['Group',sprintf('%03d',C(iC,1)),'vsGroup',sprintf('%03d',C(iC,2))],'_ORcilow.nii']);
                dvaltemp = zeros(Vraw.dim);
                dvaltemp(INDDATA) = or_twohigh(:,iC);
                DynamicBC_write_NIFTI(dvaltemp,Vraw,[VoxelDirChi2test,['Group',sprintf('%03d',C(iC,1)),'vsGroup',sprintf('%03d',C(iC,2))],'_ORcihigh.nii']);
            end
            
            %             [p, Q]= chi2test(x);
            %             [orval orcilow orcihigh]=Chi2testOR(A,B,C,D);
        end
        
        if Para.Heat.StPermutation
            
            VoxelDirPerm = [HeatDir{i},'VoxelWise_Permutation',filesep];
            VoxelDirPermcluster = [HeatDir{i},'VoxelWise_Permutationcluster',filesep];
            mkdir(VoxelDirPerm);
            mkdir(VoxelDirPermcluster);
            C = nchoosek(1:length(SetUpPara.ParaOut(i).datalab),2);
            for j = 1:size(C,1)
                IDtemp1 = SetUpPara.ParaOut(i).datalab{C(j,1)};
                IDtemp2 = SetUpPara.ParaOut(i).datalab{C(j,2)};
                DATtemp1 = full(D(:,IDtemp1));
                DATtemp2 = full(D(:,IDtemp2));
                NUM1 = length(IDtemp1);
                NUM2 = length(IDtemp2);
                [pval ClusNum] = Clinc_permtest(DATtemp1,DATtemp2,NUM1,NUM2,1,V(1),VoxelDirPermcluster);
                Pvals = reshape(pval,V(1).dim);
                DynamicBC_write_NIFTI(Pvals,V(1),[VoxelDirPerm,['Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2))],'_PermP.nii'])
                save([VoxelDirPerm,'Group',sprintf('%03d',C(j,1)),'vsGroup',sprintf('%03d',C(j,2)),'_CorrectClusNum.mat'],'ClusNum');
            end
        end
    end
    
    if Para.Heat.R
        ROIDir = [HeatDir{i},'ROIWise',filesep];
        if Para.Heat.StChisquare
            ROIDirChi2test = [HeatDir{i},'ROIWise_Chi2test',filesep];
            mkdir(ROIDirChi2test);
            for iROI = 1:size(AtlasTempU,1)
                dirROI = [ROIDir,AtlasTempU{iROI,3},filesep];
                %
                ROIDirChi2testTemp = [ROIDirChi2test,AtlasTempU{iROI,3},filesep];
                mkdir(ROIDirChi2testTemp);
                %
                markedAtlas = dir([dirROI,'Group*-markedAtlas.csv']);
                markedBoth = dir([dirROI,'Group*-markedBoth.csv']);
                markedROI = dir([dirROI,'Group*-markedROI.csv']);
                clear MAtlas MBoth MROI Nums
                for igroup = 1:length(markedAtlas)
                    aa = markedAtlas(igroup).name;
                    indsep = find(aa=='-');
                    Nums(igroup,1) = str2num(aa(indsep(end-1)+6:indsep(end)-1));
                    MAtlas(igroup,:) = load([dirROI,markedAtlas(igroup).name]);
                    MBoth(igroup,:) = load([dirROI,markedBoth(igroup).name]);
                    MROI(igroup,:) = load([dirROI,markedROI(igroup).name]);
                end
                C = nchoosek(1:length(markedAtlas),2);
                
                MAtlas2 = repmat(Nums,1,size(MAtlas,2))-MAtlas;
                MBoth2 = repmat(Nums,1,size(MBoth,2))-MBoth;
                MROI2 = repmat(Nums,1,size(MROI,2))-MROI;
                clear p_ANOVA_MAtlas Q_ANOVA_MAtlas p_two_MAtlas q_two_MAtlas or_two_MAtlas or_twolow_MAtlas or_twohigh_MAtlas
                clear p_ANOVA_MBoth Q_ANOVA_MBoth p_two_MBoth q_two_MBoth or_two_MBoth or_twolow_MBoth or_twohigh_MBoth
                clear p_ANOVA_MROI Q_ANOVA_MROI p_two_MROI q_two_MROI or_two_MROI or_twolow_MROI or_twohigh_MROI
                for idat = 1:size(MAtlas,2)
                    x = [MAtlas(:,idat),MAtlas2(:,idat)];
                    [p_ANOVA_MAtlas(idat),Q_ANOVA_MAtlas(idat),p_two_MAtlas(idat,:),q_two_MAtlas(idat,:),or_two_MAtlas(idat,:),or_twolow_MAtlas(idat,:),or_twohigh_MAtlas(idat,:)] = ClincFreq_Chi2test(x,C);
                    x = [MBoth(:,idat),MBoth2(:,idat)];
                    [p_ANOVA_MBoth(idat),Q_ANOVA_MBoth(idat),p_two_MBoth(idat,:),q_two_MBoth(idat,:),or_two_MBoth(idat,:),or_twolow_MBoth(idat,:),or_twohigh_MBoth(idat,:)] = ClincFreq_Chi2test(x,C);
                    x = [MROI(:,idat),MROI2(:,idat)];
                    [p_ANOVA_MROI(idat),Q_ANOVA_MROI(idat),p_two_MROI(idat,:),q_two_MROI(idat,:),or_two_MROI(idat,:),or_twolow_MROI(idat,:),or_twohigh_MROI(idat,:)] = ClincFreq_Chi2test(x,C);
                end
                [VROI DROI] = Dynamic_read_dir_NIFTI(AtlasTempU{iROI,1});
                
                if size(C,1)>1
                    csvwrite([ROIDirChi2testTemp,'p_ANOVA_markedAtlas.csv'],p_ANOVA_MAtlas);
                    csvwrite([ROIDirChi2testTemp,'Q_ANOVA_markedAtlas.csv'],Q_ANOVA_MAtlas);
                    csvwrite([ROIDirChi2testTemp,'p_ANOVA_markedBoth.csv'],p_ANOVA_MBoth);
                    csvwrite([ROIDirChi2testTemp,'Q_ANOVA_markedBoth.csv'],Q_ANOVA_MBoth);
                    csvwrite([ROIDirChi2testTemp,'p_ANOVA_markedROI.csv'],p_ANOVA_MROI);
                    csvwrite([ROIDirChi2testTemp,'Q_ANOVA_markedROI.csv'],Q_ANOVA_MROI);
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = p_ANOVA_MAtlas(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'p_ANOVA_markedAtlas.nii']);
                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Q_ANOVA_MAtlas(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'Q_ANOVA_markedAtlas.nii']);
                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = p_ANOVA_MBoth(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'p_ANOVA_markedBoth.nii']);
                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Q_ANOVA_MBoth(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'Q_ANOVA_markedBoth.nii']);
                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = p_ANOVA_MROI(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'p_ANOVA_markedROI.nii']);
                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Q_ANOVA_MROI(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'q_ANOVA_markedROI.nii']);
                end
                for iroi = 1:size(C,1)
                    %
                    csvwrite([ROIDirChi2testTemp,'p_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.csv'],p_two_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = p_two_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'p_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.nii']);
                    csvwrite([ROIDirChi2testTemp,'Q_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.csv'],q_two_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = q_two_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'Q_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.nii']);
                    %
                    csvwrite([ROIDirChi2testTemp,'p_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.csv'],p_two_MBoth(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = p_two_MBoth(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'p_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.nii']);
                    csvwrite([ROIDirChi2testTemp,'Q_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.csv'],q_two_MBoth(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = q_two_MBoth(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'Q_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.nii']);
                    %
                    csvwrite([ROIDirChi2testTemp,'p_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.csv'],p_two_MROI(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = p_two_MROI(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'p_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.nii']);
                    csvwrite([ROIDirChi2testTemp,'Q_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.csv'],q_two_MROI(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = q_two_MROI(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'Q_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.nii']);
                    
                    
                    csvwrite([ROIDirChi2testTemp,'OR_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.csv'],or_two_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_two_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'OR_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.nii']);
                    csvwrite([ROIDirChi2testTemp,'OR_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.csv'],or_two_MBoth(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_two_MBoth(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'OR_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.nii']);
                    csvwrite([ROIDirChi2testTemp,'OR_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.csv'],or_two_MROI(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_two_MROI(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'OR_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.nii']);
                    
                    csvwrite([ROIDirChi2testTemp,'ORcilow_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.csv'],or_twolow_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_twolow_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'ORcilow_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.nii']);
                    csvwrite([ROIDirChi2testTemp,'ORcilow_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.csv'],or_twolow_MBoth(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_twolow_MBoth(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'ORcilow_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.nii']);
                    csvwrite([ROIDirChi2testTemp,'ORcilow_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.csv'],or_twolow_MROI(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_twolow_MROI(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'ORcilow_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.nii']);
                    
                    csvwrite([ROIDirChi2testTemp,'ORcihigh_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.csv'],or_twohigh_MAtlas(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_twohigh_MAtlas(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'ORcihigh_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedAtlas.nii']);
                    csvwrite([ROIDirChi2testTemp,'ORcihigh_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.csv'],or_twohigh_MBoth(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_twohigh_MBoth(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'ORcihigh_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedBoth.nii']);
                    csvwrite([ROIDirChi2testTemp,'ORcihigh_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.csv'],or_twohigh_MROI(:,iroi));
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = or_twohigh_MROI(idats,iroi);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirChi2testTemp,'ORcihigh_Group',sprintf('%03d',C(iroi,1)),'vsGroup',sprintf('%03d',C(iroi,2)),'_markedROI.nii']);
                end
            end
        end
        
        if Para.Heat.StPermutation
            ROIDirPerm = [HeatDir{i},'ROIWise_Permutation',filesep];
            mkdir(ROIDirPerm);
            for iROI = 1:size(AtlasTempU,1)
                ROIdirss = [ROIDir,AtlasTempU{iROI,3},filesep];
                mats = dir([ROIdirss,'Single*.mat']);
                ROIDirPermU = [ROIDirPerm,AtlasTempU{iROI,3},filesep];
                mkdir(ROIDirPermU);
                
                [VROI DROI] = Dynamic_read_dir_NIFTI(AtlasTempU{iROI,1});
                
                C = nchoosek(1:length(mats),2);
                for iC = 1:size(C,1)
                    G1dat = load([ROIdirss,mats(C(iC,1)).name]);
                    G2dat = load([ROIdirss,mats(C(iC,2)).name]);
                    DATTemp1 = G1dat.IND_involve1;
                    DATTemp2 = G2dat.IND_involve1;
                    NUM1 = size(DATTemp1,1);
                    NUM2 = size(DATTemp2,1);
                    Permpval_MAtlas = Clinc_permtest(DATTemp1',DATTemp2',NUM1,NUM2,2,[],[]);
                    DATTemp1 = G1dat.IND_involve2;
                    DATTemp2 = G2dat.IND_involve2;
                    NUM1 = size(DATTemp1,1);
                    NUM2 = size(DATTemp2,1);
                    Permpval_MBoth = Clinc_permtest(DATTemp1',DATTemp2',NUM1,NUM2,2,[],[]);
                    DATTemp1 = G1dat.IND_involve3;
                    DATTemp2 = G2dat.IND_involve3;
                    NUM1 = size(DATTemp1,1);
                    NUM2 = size(DATTemp2,1);
                    Permpval_MROI = Clinc_permtest(DATTemp1',DATTemp2',NUM1,NUM2,2,[],[]);
                    
                    csvwrite([ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedAtlas.csv'],Permpval_MAtlas);                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Permpval_MAtlas(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedAtlas.nii']);
                    
                    csvwrite([ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedBoth.csv'],Permpval_MBoth);                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Permpval_MBoth(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedBoth.nii']);
                    
                    csvwrite([ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedROI.csv'],Permpval_MROI);                    
                    DAT = zeros(VROI.dim);
                    indtemp = AtlasTempU{iROI,6};
                    for idats = 1:length(indtemp)
                        DAT(indtemp{idats}) = Permpval_MROI(idats);
                    end
                    DynamicBC_write_NIFTI(DAT,VROI,[ROIDirPermU,'Group_',sprintf('%03d',C(iC,1)),'vsGroup_',sprintf('%03d',C(iC,2)),'_markedROI.nii']);
                    clear Permpval_MAtlas Permpval_MBoth Permpval_MROI
                end
            end
        end
    end
end

end