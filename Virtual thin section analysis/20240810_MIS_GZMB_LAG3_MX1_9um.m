bfCheckJavaPath;
% LAG3
metadata = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\LAG3spots.tif');
LAG3= extractChannels(1,metadata,0,[]);
LAG3_9um = LAG3(:,:,73:105);

% MX1
metadata = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\MIS_MX1_spotmask.tif');
MX1= extractChannels(1,metadata,0,[]);
MX1_9um = MX1(:,:,73:105);

% GZMB
metadata = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\MIS_GZMB_spotmask2.tif');
GZMB= extractChannels(1,metadata,0,[]);
GZMB_9um = GZMB(:,:,73:105);

% metadata = bfGetReader('/n/scratch/users/y/yil418/F8ll/Dataset1-LSP13626-melanoma_in-situ_mask.tif');
metadata = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\aggregated_full_var1_label_lower-thresh_flowremove_postrefine_resize.tif');
mask = extractChannels(1,metadata,0,[]);
mask_9um = mask(:,:,73:105);
mask_resize = imresize3(mask_9um,size(GZMB_9um),'nearest');

stats = regionprops3(mask_resize, GZMB_9um,'MeanIntensity','Volume');
writetable(stats,'N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\F8ll_results\In_situ_9_micron\Dataset1-LSP13626-melanoma_in-situ_GZMB_73_105_9um_stats.csv');

stats = regionprops3(mask_resize, LAG3_9um,'MeanIntensity','Volume');
writetable(stats,'N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\F8ll_results\In_situ_9_micron\Dataset1-LSP13626-melanoma_in-situ_LAG3_73_105_9um_stats.csv');

stats = regionprops3(mask_resize, MX1_9um,'MeanIntensity','Volume');
writetable(stats,'N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\F8ll_results\In_situ_9_micron\Dataset1-LSP13626-melanoma_in-situ_MX1_73_105_9um_stats.csv');