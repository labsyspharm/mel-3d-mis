mask = volumeRead(['F:\F8iisegmentation\aggregated_full_var1_label_lower-thresh_flowremove_postrefine_resize.tif']);
metadata =bfGetReader(['D:\F8IIaBG.tif']);
% pheno = readtable('G:\F8IIregistered\output\phenotypes.csv');
% mask = imresize3(mask,[metadata.getSizeY metadata.getSizeX metadata.getSizeZ],'nearest');
mask = imresize3(mask(509:1522,584:1715,:),size(I),'nearest');
% 
vol = [];
stats = [];


%%%
channels = [2 4 7 8 10 14 15 19 20 22 23 26 28 30 31 32 35 36 38 39 42 43 47 49 54 57 61];
channels = [11];
bgSub = [0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];
allStats=[];
allCells = [];

mask = extractChannels(1,bfGetReader(['F:\F8iisegmentation\aggregated_full_var1_label_lower-thresh_flowremove_postrefine_resize.tif']),0,[]);

for iChan=channels
    I = extractChannels(iChan,metadata,0,[]);
%     if ~(iChan == 14 || iChan == 15 || iChan == 23 || iChan == 26)
%         I=I-40;
%     end
    Icrop  =  imresize3(I,size(mask),'nearest');
    clear I
    stats= (regionprops3(mask,Icrop,'MeanIntensity','VoxelIdxList'));
    allStats= cat(2,allStats,stats.MeanIntensity);
    cells = zeros(size(mask),'uint16');
    
    for i = 1:size(stats,1)
        if ~isnan(stats.MeanIntensity(i))
            cells(stats.VoxelIdxList{i}) = stats.MeanIntensity(i)*5;
            disp(int2str(i))
        end
    end
    
    allCells=cat(5,allCells,imresize(cells,0.5,'nearest'),imresize(Icrop,0.5));
    disp(['Finished channel ' num2str(iChan)])
end
% imarisShowArray(allCells)

%% CD103 counting in MIS
metadata =bfGetReader(['D:\F8IIaBG.tif']);
mask = extractChannels(1,bfGetReader(['F:\F8iisegmentation\aggregated_full_var1_label_lower-thresh_flowremove_postrefine_resize.tif']),0,[]);
channels = [4 14 35 49];
allStats=[];
allCells=[];
for iChan=channels
    I = extractChannels(iChan,metadata,0,[]);
    Icrop  =  imresize3(I,size(mask));
    clear I
    stats= (regionprops3(mask,Icrop,'MeanIntensity','VoxelIdxList'));
    allStats= cat(2,allStats,stats.MeanIntensity);
    allCells=cat(5,allCells,imresize(Icrop,0.5));
    disp(['Finished channel ' num2str(iChan)])
end

test= imgaussfilt3(allCells(:,:,:,:,1),2);
epimask = test > thresholdOtsu(test); 
test= imgaussfilt3(allCells(:,:,:,:,2),2);
epimask=epimask+(test>thresholdOtsu(test));
epimask = imdilate(bwareaopen(epimask>0,10000),strel('disk',5));

stats= (regionprops3(mask,imresize3(epimask,size(mask)),'MeanIntensity'));

quantMIS=readtable('D:\F8iia-quantification4.csv');
epidermis_cd103 = numel(find([quantMIS.CD103>200] & [stats.MeanIntensity>0]));
dermis_cd103 = numel(find([quantMIS.CD103>200] & [stats.MeanIntensity==0]));
cd103 = numel(find([quantMIS.CD103>200]));
epidermis_cd103/cd103
dermis_cd103/cd103

%% nuclear actin fibers
metadata = loci.formats.Memoizer(bfGetReader());
metadata.setId('D:\F8IIaBG.tif') 
channels = [3, 4];
allCells = [];
for iChan=channels
    I = extractChannels(iChan,metadata,0,[]);
    Icrop  =  imresize3(I,size(mask),'nearest');
    clear I
    allCells=cat(5,allCells,Icrop);
    disp(['Finished channel ' num2str(iChan)])
end





pheno=readtable('C:\Users\cy101\Dropbox (HMS)\2023_3D (1)\data\phenotypes_MIS.csv');
quantMIS=readtable('D:\F8iia-quantification4.csv');
actin = [4317 4318 2518 2256 1978 3258 1329 2722 2313 10293 3586 1049 8243 2774 2634 7970 6183 8784 9172 8653 8243 5732 5230 6223 8005 8848 9935 9949 9351 2082 809 605 4641 8844 3803 4774 5586 919 764 6386 3369 6081 7021 6391 9531 9528 8460 5257 7727 5603 4005 4004 8042 9014 4013 4660 444 6410 10063 10453 7899 8634 9439 9590 6188 6652 8270 8678 11060 11237 10322 9065 2054 2336 2605 2599 3164 2316 10904 10298 7475 7772 7762 7332 5965 4965 4433 5040 4117 3198 3920 1658 4395 1867 1286 2527 1519 629 2111 6233 910 1261 4761 4898 10729 10329 2481 7079 2884 6755 8222 2868 1223 5636 2562 2434 1617 1012 3669 7154 1084 1961 1838];
stats= (regionprops3(mask,'VoxelIdxList','Centroid','Volume'));
cells = zeros(size(mask),'uint8');
for i = actin
    cells(stats.VoxelIdxList{i}) = 255;
    disp(int2str(i))
end

summary=[];
nearbyCD4=zeros(size(mask),'uint8');
nearbyCD8=zeros(size(mask),'uint8');
nearbydendritic=zeros(size(mask),'uint8');
nearbymacrophage=zeros(size(mask),'uint8');
index=[];
for iCell = actin
    test = mask((stats.Centroid(iCell,2)-80):(stats.Centroid(iCell,2)+80),...
                (stats.Centroid(iCell,1)-80):(stats.Centroid(iCell,1)+80),:);
    Idist = bwdist(test==iCell);
    nearbyCells=zeros(size(test));
    diststats= (regionprops3(test,Idist,'MinIntensity','MaxIntensity','VoxelIdxList'));
    CD4=0; CD4GZMB=0; CD8a=0; CD8aGZMB=0; dendritic=0; dendriticGZMB=0; totalCells=0; macrophage=0; TissueTGZMB=0; TissueT=0;
    for nearbyCell = 1:size(diststats,1)
        if isempty(diststats.MinIntensity{nearbyCell}) || nearbyCell > pheno.CellID(end,1)+1
            continue
        else
            if diststats.MinIntensity{nearbyCell}<5 %diststats.MinIntensity{nearbyCell}>25 && diststats.MinIntensity{nearbyCell}<50
                if nearbyCell == iCell
                    continue
                end
                disp('Searching')
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'CD4 T') && quantMIS.GZMB_SPOTS(nearbyCell)>0 CD4GZMB=CD4GZMB+1; end
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'CD4 T')  CD4=CD4+1; nearbyCD4(stats.VoxelIdxList{nearbyCell}) = quantMIS.CD4(nearbyCell)/max(quantMIS.CD4)*255; end
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'CD8 T') && quantMIS.GZMB_SPOTS(nearbyCell)>0  CD8aGZMB=CD8aGZMB+1; end
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'CD8 T') CD8a=CD8a+1;    end
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'Tissue T') && quantMIS.GZMB_SPOTS(nearbyCell)>0  TissueTGZMB=TissueTGZMB+1; end
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'Tissue T') TissueT=TissueT+1;  nearbyCD8(stats.VoxelIdxList{nearbyCell}) = quantMIS.CD8a(nearbyCell)/max(quantMIS.CD8a)*255; end
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'Dendritic cells') && quantMIS.GZMB_SPOTS(nearbyCell)>0  dendriticGZMB=dendriticGZMB+1; end
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'Dendritic cells') dendritic=dendritic+1; nearbydendritic(stats.VoxelIdxList{nearbyCell}) = quantMIS.CD11c(nearbyCell)/max(quantMIS.CD11c)*255; end
                if strcmp(pheno.phenotype(find(pheno.CellID==nearbyCell-1)),'Macrophage') macrophage=macrophage+1; nearbymacrophage(stats.VoxelIdxList{nearbyCell}) = quantMIS.CD206(nearbyCell)/max(quantMIS.CD206)*255; end
                totalCells = totalCells + 1;
            end
        end
    end
    summary=cat(1,summary,cat(2,CD4, CD4GZMB, CD8a, CD8aGZMB, TissueT, TissueTGZMB, dendritic, dendriticGZMB, macrophage, totalCells));
end
test1=setdiff(unique(mask(:)),actin);
[h,p]=ttest2(stats.Volume(actin),stats.Volume(test1(2:end)))

%sox10 sox9 prame
PRAMEI = extractChannels(30,metadata,0,[]);
MART1 = extractChannels(4,metadata,0,[]);
statsPRAME = regionprops3(mask,imresize3(PRAMEI,size(mask)),'MeanIntensity');
SOX10=zeros(size(mask),'uint8');
SOX9=zeros(size(mask),'uint8');
SOX109=zeros(size(mask),'uint8');
PRAME=zeros(size(mask),'uint8');
keratinocytes=zeros(size(mask),'uint8');
for iCell = 1:size(stats,1)
    if quantMIS.SOX10(iCell)> 400 && quantMIS.SOX9ForGating(iCell)<1 && quantMIS.MART1(iCell)>1200 SOX10(stats.VoxelIdxList{iCell}) = quantMIS.SOX10(iCell)/max(quantMIS.SOX10)*255; %SOX10+ SOX9-
    elseif quantMIS.SOX10(iCell)< 400 && quantMIS.SOX9ForGating(iCell)>0 && quantMIS.MART1(iCell)>1200 && quantMIS.panCK(iCell)<300 SOX9(stats.VoxelIdxList{iCell}) = quantMIS.SOX9(iCell)/max(quantMIS.SOX9)*255; %SOX10- SOX9+
    elseif quantMIS.SOX10(iCell)> 400 && quantMIS.SOX9ForGating(iCell)>0 && quantMIS.MART1(iCell)>1200 && quantMIS.panCK(iCell)<300 SOX109(stats.VoxelIdxList{iCell}) = quantMIS.SOX10(iCell)/max(quantMIS.SOX10)*255; %SOX10+ SOX9+
    end
    if statsPRAME.MeanIntensity(iCell) >250 PRAME(stats.VoxelIdxList{iCell}) = statsPRAME.MeanIntensity(iCell)/max(statsPRAME.MeanIntensity)*255; 
    end
    if quantMIS.panCK(iCell)>300 keratinocytes(stats.VoxelIdxList{iCell}) = 255;
    end

   disp(int2str(iCell))
end
% tiffwriteimj(cat(4,allCells(:,:,:,:,1),allCells(:,:,:,:,2),SOX10,SOX9,SOX109,PRAME,cells,nearbyCD4,nearbyCD8,nearbydendritic,nearbymacrophage), 'D:\MIS_actinmelanocytes_immunecells_2.tif')
tiffwriteimj(cat(4,SOX10,SOX9,SOX109,PRAME,bwperim(imresize3(keratinocytes,size(PRAME),'nearest'))),'D:\SOX10SOX9PRAME_maskplot_MIS.tif')


%% visually gate
img=cat(3,max(allCells(:,:,:,:,4),[],3)*10,max(allCells(:,:,:,:,6),[],3)*5,allCells(:,:,10,:,1));
outlines = bwperim(max(allCells(:,:,:,:,3)>500,[],3));
img(:,:,3)=outlines*65535;
imshow(img,[])
figure,imshowpair(bwperim(max(allCells(:,:,:,:,3)>500,[],3)),max(allCells(:,:,:,:,4),[],3)*10)

%%% SOX9
channels = [47];
thresh = [195];
sigma = [4];
for iChan = 1:numel(channels)
    I = extractChannels(channels(iChan),metadata,0,[]);
    Icrop  =  imresize3(I,0.25,'nearest');
    clear I
%     allCells= cat(5,allCells,Icrop);
    Ithcrop = filterLoGND(Icrop,sigma(iChan));
    Ithcrop = -(Ithcrop-max(Ithcrop(:)));
    Imax = imregionalmax(Ithcrop);
    Imin = Imax.*(Ithcrop>thresh(iChan));
%     allCells=cat(5,allCells,Imin);

    %SOX9 spot gating
    stats= (regionprops3(mask,imresize3(Imin,size(mask)),'MeanIntensity','Volume','VoxelIdxList'));
    cells = zeros(size(mask),'uint16');
    for i = 1:size(stats,1)
        if ~isnan(stats.MeanIntensity(i))
            cells(stats.VoxelIdxList{i}) = stats.MeanIntensity(i)*stats.Volume(i)*10;
            disp(int2str(i))
        end
    end
    allStats= cat(2,allStats,(stats.MeanIntensity.*stats.Volume));

%     %SOX9 mean intensity
%     stats= (regionprops3(mask,Icrop,'MeanIntensity','Volume','VoxelIdxList'));
%     cells = zeros(size(mask),'uint16');
%     for i = 1:size(stats,1)
%         if ~isnan(stats.MeanIntensity(i))
%             cells(stats.VoxelIdxList{i}) = stats.MeanIntensity(i)*10;
%             disp(int2str(i))
%         end
%     end
%     disp(['Finished channel ' num2str(iChan)])
%     allStats= cat(2,allStats,(stats.MeanIntensity));

end




% % % % collagen
% I = extractChannels(71,metadata,0);
% Icrop= imresize3(I,0.5);
% collagenmask = bwareaopen(Icrop > thresholdOtsu(Icrop),20);
% Idist = bwdist(collagenmask);
% stats=table2array(regionprops3(imresize3(mask,0.5,'nearest'),Idist,'MeanIntensity'));
% allStats= cat(2,allStats,stats);


% metadata =bfGetReader(['F:\F8iisegmentation\F8iiacollagenDistancemap.tif']);
% I = extractChannels(1,metadata,0);
% stats= cat(2,stats,table2array(regionprops3(mask,I,'MeanIntensity')));

% % centroid, volume
morphstats= table2array(regionprops3(mask,'Centroid','Volume'));
morphstats(:,2) = morphstats(:,2)*0.14;
morphstats(:,3) = morphstats(:,3)*0.14;
morphstats(:,4) = morphstats(:,4)*0.28;
allStats= cat(2,allStats,morphstats);
morphstats= table2array(regionprops3(mask,'PrincipalAxisLength'));
allmorphstats=zeros(numel(morphstats),3);
for iCell = 1:numel(morphstats)
    if morphstats{iCell,1} > 0
        allmorphstats(iCell,:) =morphstats{iCell}';
    end
end
allStats= cat(2,allStats,allmorphstats);
%X_centroid Y_centroid Z_centroid volume


% MIS spot detection
channels = [3 61 27 66];
thresh = [-220 -270 -350 -100];
sigma = [ 2 2 2 2];
Imin=[];
metadata =bfGetReader(['D:\F8IIaBG.tif']);
for iChan = 3:4%1:numel(channels)
    I = extractChannels(channels(iChan),metadata,0,[]);
    for i=1:(size(I,2)/1000+1)
        if i >10
            Icrop  =  I(:,1000*(i-1)+1:end,:);
        else
            Icrop  =  I(:,1000*(i-1)+1:1000*(i),:);
        end
        
    Ithcrop = filterLoGND(Icrop,sigma(iChan));
    Imax = imregionalmin(Ithcrop);
    Imin = cat(2,Imin,uint8(Imax.*(Ithcrop<thresh(iChan))));
    
    clear Imax Ithcrop 
    disp(['Running tile ' num2str(i)])
    end
 clear I
    stats= (regionprops3(imresize3(mask,size(Imin),'nearest'),Imin,'MeanIntensity','Volume'));
%     cells = zeros(size(maskCrop),'uint8');
%     for i = 1:size(stats,1)
%         if ~isnan(stats.MeanIntensity(i))
%             cells(stats.VoxelIdxList{i}) = stats.Volume(i)*stats.MeanIntensity(i)*10;
%             disp(int2str(i))
%         end
%     end
%     allCells=cat(5,allCells,cells);
    allStats= cat(2,allStats,(stats.MeanIntensity).*(stats.Volume));
    disp(['Finished channel ' num2str(iChan)])
end

mhc1=find(cat(1,mhc1stats.MeanIntensity>1000));
scatter(allStats(mhc1),mhc1stats.MeanIntensity(mhc1))
cd4=find(cat(1,cd4stats.MeanIntensity>300));
noncd4=find(cat(1,cd4stats.MeanIntensity<300));
sum(allStats(cd4))
cd8=find(cat(1,cd8stats.MeanIntensity>600));
noncd8=find(cat(1,cd8stats.MeanIntensity<600));
sum(allStats(noncd8))

mart1=find([mart1stats.MeanIntensity>4000] & [DEJmaskstats.MeanIntensity>0]);
scatter(allStats(mart1),mhc1stats.MeanIntensity(mart1)-100)

Icrop = I(4500:5501,3500:4501,:);
stats= (regionprops3(imresize(mask(2250:2750,750:1250,:),2,'nearest'),Imin,'MeanIntensity','Volume'));


%% IM spot detection
channels = [66];
thresh = [-100];
sigma = [ 2];
Imin=[];
metadata =bfGetReader(['V:\cycif-techdev\ClarenceLSM980\F8IIc16bitbg.tif']);
for iChan = 1:numel(channels)
    I = extractChannels(channels(iChan),metadata,0,[5000 5000 1000 1000]);
    for i=1:(size(I,2)/1000+1)
        if i >10
            Icrop  =  I(:,1000*(i-1)+1:end,:);
        else
            Icrop  =  I(:,1000*(i-1)+1:1000*(i),:);
        end
        
    Ithcrop = filterLoGND(Icrop,sigma(iChan));
    Imax = imregionalmin(Ithcrop);
    Imin = cat(2,Imin,uint8(Imax.*(Ithcrop<thresh(iChan))));
    
    clear Imax Ithcrop 
    disp(['Running tile ' num2str(i)])
    end
 clear I
    stats= (regionprops3(imresize3(mask,size(Imin),'nearest'),Imin,'MeanIntensity','Volume'));
%     cells = zeros(size(maskCrop),'uint8');
%     for i = 1:size(stats,1)
%         if ~isnan(stats.MeanIntensity(i))
%             cells(stats.VoxelIdxList{i}) = stats.Volume(i)*stats.MeanIntensity(i)*10;
%             disp(int2str(i))
%         end
%     end
%     allCells=cat(5,allCells,cells);
    allStats= cat(2,allStats,(stats.MeanIntensity).*(stats.Volume));
    disp(['Finished channel ' num2str(iChan)])
end
%% MX1 spots associating with MHC1
metadata =bfGetReader(['D:\F8IIaBG.tif']);
mhc1 = extractChannels(7,metadata,0,[]);

MX1 = extractChannels(1,bfGetReader('D:\MX1spots.tif'),0,[]);
DEJ=imread('D:\DEJmask.tif');
DEJ3D = padarray(DEJ>0,[0 0 193],'replicate','post');
MX1 = uint8(MX1).*uint8(DEJ3D);

spots = find(MX1>0);
[row,col,page] = ind2sub(size(MX1),spots)


spotsds = zeros(round(size(MX1)/4),'uint8');
for iSpot = 1:numel(row)
    spotsds(round(row(iSpot)/4),round(col(iSpot)/4),round(page(iSpot)/4))=1;
end
Idist =bwdist(spotsds);
Imin = imimposemin(Idist,spotsds);
spotmask = Idist<5;
spotsLabel = watershed(Imin);
spotsLabel = spotsLabel.*cast(spotmask,class(spotsLabel));
stats = regionprops3(spotsLabel,imresize3(mhc1,size(spotsLabel)),'MeanIntensity');
mean(stats.MeanIntensity)
std(stats.MeanIntensity)/sqrt(numel(stats))

unspotmask = imresize3(DEJ3D,size(spotmask),'nearest') & ~spotmask;
mhc1ds = imresize3(mhc1,size(unspotmask));
mean(mhc1ds(unspotmask))



x = 1:2;
data = [mean(stats.MeanIntensity) mean(mhc1ds(unspotmask))]';
errhigh = [std(stats.MeanIntensity)/sqrt(numel(stats)) 0];
errlow  = [std(stats.MeanIntensity)/sqrt(numel(stats)) 0];

bar(x,data)                

hold on

er = errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  


%% MIS Ki67
quantMIS=readtable('D:\F8iia-quantification4.csv');
CD3E_ki67 = numel(find([quantMIS.CD3E>250] & [quantMIS.Ki67>160]))
CD31_ki67 = numel(find([quantMIS.CD31>578] & [quantMIS.Ki67>160]))
monocyte_ki67= numel(find([quantMIS.CD11c>120] & [quantMIS.CD11b<300] & [quantMIS.Ki67>160]))
Ki67=  numel(find([quantMIS.Ki67>160]))
CD3E_ki67/Ki67*100
CD31_ki67/Ki67*100
monocyte_ki67/Ki67*100

%% IM Ki67
quantIM=readtable('C:\Users\cy101\Dropbox (HMS)\2023_3D (1)\data\F8iic-quantificationV2.csv');
CD3E_ki67 = numel(find([quantIM.CD3E>240] & [quantIM.Ki67>100]))
tumor_ki67 = numel(find([quantIM.MART1>800] & [quantIM.Ki67>100]))
Ki67=  numel(find([quantIM.Ki67>100]))
monocyte_ki67= numel(find([quantIM.CD11c>100] & [quantIM.CD11b<220] & [quantIM.Ki67>100]))
macrophage_ki67= numel(find([quantIM.CD163>700] & [quantIM.Ki67>100]))
CD8_ki67= numel(find([quantIM.CD8a>100] & [quantIM.Ki67>100]))


tumor_ki67/Ki67*100
monocyte_ki67/Ki67*100
CD3E_ki67/Ki67*100
macrophage_ki67/Ki67*100

%% IM GZMB


quantIM=readtable('C:\Users\cy101\Dropbox (HMS)\2023_3D (1)\data\F8iic-quantificationV2.csv');
quantMIS=readtable('D:\F8iia-quantification4.csv');

CD4_GZMB_MIS = histogram(quantMIS.GZMB_SPOTS((find([quantMIS.CD4>1050] & [quantMIS.GZMB_SPOTS>2]))));
numel(find([quantMIS.CD4>1050] & [quantMIS.GZMB_SPOTS>2]))
mean(quantMIS.GZMB_SPOTS((find([quantMIS.CD4>1050] & [quantMIS.GZMB_SPOTS>2]))))
hold on
CD8_GZMB_MIS = histogram(quantMIS.GZMB_SPOTS((find([quantMIS.CD8a>450] & [quantMIS.GZMB_SPOTS>2]))));
numel(find([quantMIS.CD8a>450] & [quantMIS.GZMB_SPOTS>2]))
mean(quantMIS.GZMB_SPOTS((find([quantMIS.CD8a>450] & [quantMIS.GZMB_SPOTS>2]))))

CD4_GZMB_IM = histogram(quantIM.GZMB_spots(find([quantIM.CD4>800] & [quantIM.GZMB_spots>2])));
numel(find([quantIM.CD4>800] & [quantIM.GZMB_spots>2]))
mean(quantIM.GZMB_spots(find([quantIM.CD4>800] & [quantIM.GZMB_spots>2])))
hold on
CD8_GZMB_IM = histogram(quantIM.GZMB_spots(find([quantIM.CD8a>400] & [quantIM.GZMB_spots>2])));
numel(find([quantIM.CD8a>400] & [quantIM.GZMB_spots>2]))
mean(quantIM.GZMB_spots(find([quantIM.CD8a>400] & [quantIM.GZMB_spots>2])))

CD20_IM = extractChannels(28,metadata,0,[]);
CD4_IM = extractChannels(26,metadata,0,[]);
CD8_IM = extractChannels(36,metadata,0,[]);
CD11c_IM = extractChannels(43,metadata,0,[]);
MX1_IM = extractChannels(3,metadata,0,[]);
GZMB_IM = extractChannels(66,metadata,0,[]);
collagen_IM = extractChannels(70,metadata,0,[]);
CD20 = [11010 8410 7924 8744 1299 3241 3517 6398 4164 6536 5941 2084 2946 551 3686 7158]  ;
stats= (regionprops3(mask,'Centroid'));
allCells = zeros(400,400,size(CD20_IM,3),numel(CD20),7,'uint8');
for iCell = 1:numel(CD20)

   allCells(:,:,:,iCell,1) = CD20_IM(stats.Centroid(CD20(iCell),2)*2-199:stats.Centroid(CD20(iCell),2)*2+200,...
                            stats.Centroid(CD20(iCell),1)*2-199:stats.Centroid(CD20(iCell),1)*2+200,:);
   allCells(:,:,:,iCell,2) = CD4_IM(stats.Centroid(CD20(iCell),2)*2-199:stats.Centroid(CD20(iCell),2)*2+200,...
                            stats.Centroid(CD20(iCell),1)*2-199:stats.Centroid(CD20(iCell),1)*2+200,:);
   allCells(:,:,:,iCell,3) = CD8_IM(stats.Centroid(CD20(iCell),2)*2-199:stats.Centroid(CD20(iCell),2)*2+200,...
                            stats.Centroid(CD20(iCell),1)*2-199:stats.Centroid(CD20(iCell),1)*2+200,:);
   allCells(:,:,:,iCell,4) = CD11c_IM(stats.Centroid(CD20(iCell),2)*2-199:stats.Centroid(CD20(iCell),2)*2+200,...
                            stats.Centroid(CD20(iCell),1)*2-199:stats.Centroid(CD20(iCell),1)*2+200,:);
   allCells(:,:,:,iCell,5) = MX1_IM(stats.Centroid(CD20(iCell),2)*2-199:stats.Centroid(CD20(iCell),2)*2+200,...
                            stats.Centroid(CD20(iCell),1)*2-199:stats.Centroid(CD20(iCell),1)*2+200,:);
    allCells(:,:,:,iCell,6) = GZMB_IM(stats.Centroid(CD20(iCell),2)*2-199:stats.Centroid(CD20(iCell),2)*2+200,...
                            stats.Centroid(CD20(iCell),1)*2-199:stats.Centroid(CD20(iCell),1)*2+200,:);
   allCells(:,:,:,iCell,7) = collagen_IM(stats.Centroid(CD20(iCell),2)*2-199:stats.Centroid(CD20(iCell),2)*2+200,...
                            stats.Centroid(CD20(iCell),1)*2-199:stats.Centroid(CD20(iCell),1)*2+200,:);
end
    


LAG3_IM = extractChannels(27,metadata,0,[5000 5000 1000 1000]);
GZMB_IM = extractChannels(66,metadata,0,[5000 5000 1000 1000]);
GZMBspot_IM = extractChannels(1,bfGetReader('D:\IM_GZMB_spotmask.tif'),0,[5000 5000 1000 1000]);
mask_IM = extractChannels(1,bfGetReader('D:\F8IIc IM mask.tif'),0,[]);


CD8_IM = extractChannels(36,metadata,0,[]);
CD103_IM = extractChannels(49,metadata,0,[]);
CD3E_IM = extractChannels(35,metadata,0,[]);
FOXP3_IM = extractChannels(39,metadata,0,[]);


stats= (regionprops3(mask_IM,'VoxelIdxList'));
pheno = readtable('C:\Users\cy101\Dropbox (HMS)\2023_3D (1)\data\phenotypes.csv');
phenostats=find([strcmp(pheno.phenotype,'T cells')]) ; CellIDs = pheno.CellID(phenostats)+1;

 cells = zeros(size(mask_IM),'uint8');

    for i = CellIDs'
%         if ismember(i,maskIDs)
            cells(stats.VoxelIdxList{i}) = 255;
            disp(int2str(i))
%         end
    end
