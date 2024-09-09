% Dataset1-LSP13626-invasive_margin.tif
bfCheckJavaPath;
metadata_1 = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin.tif');
PDL1= extractChannels(16,metadata_1,0,[]); 

% Mask
metadata_2 = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_mask.tif');
mask = extractChannels(1,metadata_2,0,[]);
mask_resize = imresize3(mask,size(PDL1),'nearest')
stats = regionprops3(mask_resize, PDL1,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_PDL1_stats.csv');

% PDL1
PDL1_9um = PDL1(:,:,73:105)
mask_9um = mask(:,:,73:105);
mask_resize = imresize3(mask_9um,size(PDL1_9um),'nearest');
stats = regionprops3(mask_resize, PDL1_9um,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_PDL1_73_105_9um_stats.csv');

% CD103
CD103= extractChannels(49,metadata_1,0,[]); 
CD103_9um = CD103(:,:,73:105)
stats = regionprops3(mask_resize, CD103_9um,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_CD103_73_105_9um_stats.csv');

% MART1
MART1= extractChannels(4,metadata_1,0,[]); 
MART1_9um = MART1(:,:,73:105)
stats = regionprops3(mask_resize, MART1_9um,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_MART1_73_105_9um_stats.csv');

% PD1
PD1= extractChannels(40,metadata_1,0,[]); 
PD1_9um = PD1(:,:,73:105)
stats = regionprops3(mask_resize, PD1_9um,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_PD1_73_105_9um_stats.csv');


% Dataset1-LSP13626-melanoma-in-situ.tif
bfCheckJavaPath;
metadata_1 = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ.tif');
PDL1= extractChannels(16,metadata_1,0,[]);

% Mask
metadata_2 = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_mask.tif');
mask = extractChannels(1,metadata_2,0,[]);
mask_resize = imresize3(mask,size(PDL1),'nearest');
stats = regionprops3(mask_resize, PDL1,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_PDL1_stats.csv');

% PDL1
PDL1_9um = PDL1(:,:,73:105)
mask_9um = mask(:,:,73:105);
mask_resize = imresize3(mask_9um,size(PDL1_9um),'nearest');
stats = regionprops3(mask_resize, PDL1_9um,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_PDL1_73_105_9um_stats.csv');

% CD103
CD103= extractChannels(49,metadata_1,0,[]); 
CD103_9um = CD103(:,:,73:105)
stats = regionprops3(mask_resize, CD103_9um,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_CD103_73_105_9um_stats.csv');

% MART1
MART1= extractChannels(4,metadata_1,0,[]); 
MART1_9um = MART1(:,:,73:105)
stats = regionprops3(mask_resize, MART1_9um,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MART1_73_105_9um_stats.csv');

% PD1
PD1= extractChannels(40,metadata_1,0,[]); 
PD1_9um = PD1(:,:,73:105)
stats = regionprops3(mask_resize, PD1_9um,'MeanIntensity');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_PD1_73_105_9um_stats.csv');
