bfCheckJavaPath;
metadata = bfGetReader('Z:\Dataset1-LSP13626-invasive_margin.tif');
% metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin.tif');
% Extract MHC-I channel
MHCI = extractChannels(7,metadata,0,[]);
% Extract Hoechst channel
% Hoechst = extractChannels(1,metadata,0,[]);

% Downsample if needed, resize 3D volumetric intensity inage "imresize3"
MHCI_resize = imresize3(MHCI,[4792 5440 194],'nearest');
% 3-D gausssian filtering of 3D images  "imgaussfilt3"
MHCI_blur = imgaussfilt3(MHCI_resize,0.5);
% Global image threshold using 'Otsu' method  "graythresh"
threshold = graythresh(MHCI_blur);
MHCI_binary = imbinarize(MHCI_blur,threshold);
% Remove small objects from binary image "bwareaopen"
MHCI_bw = bwareaopen(MHCI_binary,50);
% Fill image regions and holes "imfill"
% MHCI_imf = imfill(MHCI_bw,'holes');
% Morphologically close image "imclose"
se = strel('disk',10);
MHCI_imc = imclose(MHCI_bw,se);
% 	imshow(max(MHCI_resize,[],3))
% 	imshowpair(max(MHCI_blur,[],3),max(MHCI_binary,[],3),'montage');
% 	imshowpair(max(MHCI_binary,[],3),max(MHCI_bw,[],3),'montage');
% 	imshowpair(max(MHCI_imf,[],3),max(MHCI_imc,[],3),'montage');

% Distance transform of a binary image "bwdist"
MHCI_distransform = bwdist(MHCI_imc);
% imshowpair(MHCI_bw(:,:,25),MHCI_distransform(:,:,25),'montage');
% imshowpair(max(MHCI_bw,[],3),max(MHCI_distransform,[],3),'montage');

% Add image label mask
metadata = bfGetReader('Z:\segmentation\F8IIc_IM_mask.tif');
% metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_mask.tif');
Mask = extractChannels(1,metadata,0,[]);

% Measure properties of 3-D volumetric image regions "regionprops3"
stats = regionprops3(Mask,MHCI_distransform,'MeanIntensity','Volume','Centroid');
writetable(stats,'Z:\F8ll_results\VGP_distance_transform.csv');