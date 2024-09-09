% LAG3
bfCheckJavaPath;
metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/IM_LAG3spots.tif');
LAG3= extractChannels(1,metadata,0,[]);

metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_mask.tif');
mask = extractChannels(1,metadata,0,[]);
mask_resize = imresize3(mask,size(LAG3),'nearest');

stats = regionprops3(mask_resize, LAG3,'MeanIntensity','Volume');
writetable(stats,'/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-invasive_margin_LAG3_stats.csv');
