% Dataset1-LSP13626-melanoma_in-situ.tif
bfCheckJavaPath;
metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ.tif');
MHCI= extractChannels(7,metadata,0,[]);
metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_mask.tif');
mask = extractChannels(1,metadata,0,[]);

% 18um section
MHCI_18um = MHCI(:,:,74:138); 
mask_18um = mask(:,:,74:138);
mask_resize = imresize3(mask_18um,size(MHCI_18um),'nearest');
stats = regionprops3(mask_resize, MHCI_18um,'MeanIntensity','Volume','Centroid','BoundingBox');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MHCI_74_138_18um_stats.csv');

% 27um section
MHCI_27um = MHCI(:,:,56:152); 
mask_27um = mask(:,:,56:152);
mask_resize = imresize3(mask_27um,size(MHCI_27um),'nearest');
stats = regionprops3(mask_resize, MHCI_27um,'MeanIntensity','Volume','Centroid','BoundingBox');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MHCI_56_152_27um_stats.csv');

% 36um section
MHCI_36um = MHCI(:,:,40:171); 
mask_36um = mask(:,:,40:171);
mask_resize = imresize3(mask_36um,size(MHCI_36m),'nearest');
stats = regionprops3(mask_resize, MHCI_36um,'MeanIntensity','Volume','Centroid','BoundingBox');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MHCI_40_171_36um_stats.csv');

% 9um section
MHCI_9um = MHCI(:,:,73:105);   
mask_9um = mask(:,:,73:105);
mask_resize = imresize3(mask_9um,size(MHCI_9um),'nearest');
stats = regionprops3(mask_resize, MHCI_9um,'MeanIntensity','Volume','Centroid','BoundingBox');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MHCI_73_105_9um_stats.csv');

% 9um section
MHCI_9um = MHCI(:,:,40:72);   
mask_9um = mask(:,:,40:72);
mask_resize = imresize3(mask_9um,size(MHCI_9um),'nearest');
stats = regionprops3(mask_resize, MHCI_9um,'MeanIntensity','Volume','Centroid','BoundingBox');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MHCI_40_72_9um_stats.csv');

% 9um section
MHCI_9um = MHCI(:,:,106:138);  
mask_9um = mask(:,:,106:138);
mask_resize = imresize3(mask_9um,size(MHCI_9um),'nearest');
stats = regionprops3(mask_resize, MHCI_9um,'MeanIntensity','Volume','Centroid','BoundingBox');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MHCI_106_138_9um_stats.csv');

% 9um section
MHCI_9um = MHCI(:,:,139:171);   
mask_9um = mask(:,:,139:171);
mask_resize = imresize3(mask_9um,size(MHCI_9um),'nearest');
stats = regionprops3(mask_resize, MHCI_9um,'MeanIntensity','Volume','Centroid','BoundingBox');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MHCI_139_171_9um_stats.csv');


