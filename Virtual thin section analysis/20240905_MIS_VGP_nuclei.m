% Dataset1-LSP13626-invasive_margin.tif
bfCheckJavaPath;
metadata_1 = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\F8iia_nucleiMask_1.tif');
nuclei_MIS_mask = extractChannels(1,metadata_1,0,[]);

metadata_2 = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\F8iic_nucleiMask_1.tif');
nuclei_VGP_mask = extractChannels(1,metadata_2,0,[]);

metadata_3 = bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\segmentation\F8lla_MIS_mask.tif');
cell_mask = extractChannels(1,metadata_3,0,[]);

nuclei_MIS_mask = imbinarize(nuclei_MIS_mask);
nuclei_VGP_mask = imbinarize(nuclei_VGP_mask);



% 9um section  
nuclei_MIS_mask_9um = nuclei_MIS_mask(:,:,73:105);
nuclei_VGP_mask_9um = nuclei_VGP_mask(:,:,73:105);

% 18um section
nuclei_MIS_mask_18um = nuclei_MIS_mask(:,:,74:138);
nuclei_VGP_mask_18um = nuclei_VGP_mask(:,:,74:138);

% 27um section
nuclei_MIS_mask_27um = nuclei_MIS_mask(:,:,56:152);
nuclei_VGP_mask_27um = nuclei_VGP_mask(:,:,56:152);

nuclei_MIS_mask_35um = nuclei_MIS_mask(:,:,40:171);
nuclei_VGP_mask_35um = nuclei_VGP_mask(:,:,40:171);

stats =regionprops3(nuclei_MIS_mask_9um,'Volume','Centroid','BoundingBox');
stats =regionprops3(nuclei_VGP_mask_9um,'Volume','Centroid','BoundingBox');

stats_18 =regionprops3(nuclei_MIS_mask_18um,'Volume','Centroid','BoundingBox');
stats_18 =regionprops3(nuclei_VGP_mask_18um,'Volume','Centroid','BoundingBox');

stats_27 =regionprops3(nuclei_MIS_mask_27um,'Volume','Centroid','BoundingBox');
stats_27 =regionprops3(nuclei_VGP_mask_27um,'Volume','Centroid','BoundingBox');

stats_35_MIS =regionprops3(nuclei_MIS_mask_35um,'Volume','Centroid','BoundingBox');
stats_35_VGP =regionprops3(nuclei_VGP_mask_35um,'Volume','Centroid','BoundingBox');

stats =regionprops3(nuclei_MIS_mask,'Volume','Centroid','BoundingBox');
stats =regionprops3(nuclei_VGP_mask,'Volume','Centroid','BoundingBox');

writetable(stats,'Z:\F8ll_results\In_situ_9_micron\Dataset1-LSP13626-melanoma_in-situ_nuclei_stats.csv');
writetable(stats,'Z:\F8ll_results\Invasive_margin_9_micron\Dataset1-LSP13626-melanoma-invasive_margin_nuclei_stats.csv');

writetable(stats_18,'Z:\F8ll_results\In_situ_9_micron\Dataset1-LSP13626-melanoma_in-situ_nuclei_18um_stats.csv');
writetable(stats_18,'Z:\F8ll_results\Invasive_margin_9_micron\Dataset1-LSP13626-melanoma-invasive_margin_nuclei_18um_stats.csv');

writetable(stats_27,'Z:\F8ll_results\In_situ_9_micron\Dataset1-LSP13626-melanoma_in-situ_nuclei_27um_stats.csv');
writetable(stats_27,'Z:\F8ll_results\Invasive_margin_9_micron\Dataset1-LSP13626-melanoma-invasive_margin_nuclei_27um_stats.csv');

writetable(stats_35_MIS,'Z:\F8ll_results\In_situ_9_micron\Dataset1-LSP13626-melanoma_in-situ_nuclei_35um_stats.csv');
writetable(stats_35_VGP,'Z:\F8ll_results\Invasive_margin_9_micron\Dataset1-LSP13626-melanoma-invasive_margin_nuclei_35um_stats.csv');