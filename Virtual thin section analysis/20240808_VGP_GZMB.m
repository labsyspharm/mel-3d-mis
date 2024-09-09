% GZMB
bfCheckJavaPath;
metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/MIS_GZMB_spotmask.tif');
GZMB= extractChannels(1,metadata,0,[]);

metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_mask.tif');
mask = extractChannels(1,metadata,0,[]);
mask_resize = imresize3(mask,size(GZMB),'nearest');

stats = regionprops3(mask_resize, GZMB,'MeanIntensity','Volume');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_GZMB_stats.csv');