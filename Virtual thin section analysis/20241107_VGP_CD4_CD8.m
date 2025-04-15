bfCheckJavaPath;
metadata_1 = bfGetReader('Z:\Dataset1-LSP13626-invasive_margin.tif');
metadata_2 = bfGetReader('Z:\segmentation\F8IIc_IM_mask.tif');
mask = extractChannels(1,metadata_2,0,[]);

% CD4
CD4= extractChannels(26,metadata_1,0,[]); 
CD4_9um = CD4(:,:,73:105);
mask_resize = imresize3(mask,size(CD4_9um),'nearest');
stats = regionprops3(mask_resize, CD4_9um,'MeanIntensity');
writetable(stats,'Z:\F8ll_results\Dataset1-LSP13626-invasive_margin_CD4_73_105_9um_stats.csv');

% CD8
CD8= extractChannels(36,metadata_1,0,[]); 
CD8_9um = CD8(:,:,73:105);
stats = regionprops3(mask_resize, CD8_9um,'MeanIntensity');
writetable(stats,'Z:\F8ll_results\Dataset1-LSP13626-invasive_margin_CD8_73_105_9um_stats.csv');