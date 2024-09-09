% LAG3
bfCheckJavaPath;
% metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/LAG3spots.tif');
metadata = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\LAG3spots.tif');
LAG3= extractChannels(1,metadata,0,[]);


% metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_mask.tif');
metadata = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\aggregated_full_var1_label_lower-thresh_flowremove_postrefine_resize.tif');
mask = extractChannels(1,metadata,0,[]);
mask_resize = imresize3(mask,size(LAG3),'nearest');

stats = regionprops3(mask_resize, LAG3,'MeanIntensity','Volume');
writetable(stats,'N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\F8ll_results\In_situ\Dataset1-LSP13626-melanoma_in-situ_LAG3_97_115_5um_stats.csv');