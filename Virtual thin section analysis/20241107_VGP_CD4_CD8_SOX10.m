bfCheckJavaPath;
metadata_1 = bfGetReader('Z:\Dataset1-LSP13626-invasive_margin.tif');
metadata_2 = bfGetReader('Z:\segmentation\F8IIc_IM_mask.tif');
mask = extractChannels(1,metadata_2,0,[]);

% CD4
CD4= extractChannels(26,metadata_1,0,[]); 
CD4_9um = CD4(:,:,73:105);
mask_resize = imresize3(mask,size(CD4),'nearest');
mask_9um = mask_resize(:,:,73:105);
stats = regionprops3(mask_9um, CD4_9um,'MeanIntensity');
writetable(stats,'Z:\F8ll_results\Dataset1-LSP13626-invasive_margin_CD4_73_105_9um_stats.csv');

% CD8
CD8= extractChannels(36,metadata_1,0,[]); 
CD8_9um = CD8(:,:,73:105);
stats = regionprops3(mask_9um, CD8_9um,'MeanIntensity');
writetable(stats,'Z:\F8ll_results\Dataset1-LSP13626-invasive_margin_CD8_73_105_9um_stats.csv');

% SOX10
metadata_3 = bfGetReader('Z:\segmentation\F8iic_nucleiMask_1.tif');
nuclei_mask =  extractChannels(1,metadata_3,0,[]);
nuclei_mask = double(mask) .* double(nuclei_mask > 0);
SOX10= extractChannels(8,metadata_1,0,[]); 
SOX10_9um = SOX10(:,:,73:105);

nuclei_mask_resize = imresize3(nuclei_mask,size(SOX10),'nearest');
nuclei_mask_9um =nuclei_mask_resize(:,:,73:105);
stats = regionprops3(nuclei_mask_9um, SOX10_9um,'MeanIntensity');
writetable(stats,'Z:\F8ll_results\Dataset1-LSP13626-invasive_margin_SOX10_73_105_9um_stats.csv');


