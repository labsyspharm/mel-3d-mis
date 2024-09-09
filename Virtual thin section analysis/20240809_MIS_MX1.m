% MX1
bfCheckJavaPath;
metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/MIS_MX1_spotmask.tif');
MX1= extractChannels(1,metadata,0,[]);

metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_mask.tif');
mask = extractChannels(1,metadata,0,[]);
mask_resize = imresize3(mask,size(MX1),'nearest');

stats = regionprops3(mask_resize, MX1,'MeanIntensity','Volume');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_MX1_stats.csv');