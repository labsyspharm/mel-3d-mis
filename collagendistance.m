imagePath = 'V:\cycif-techdev\ClarenceLSM980\20micronF8II\Dataset1-LSP13626-melanoma_in-situ.tif';
metadata = loci.formats.Memoizer(bfGetReader());
metadata.setId(imagePath)

collagen = extractChannels(54,metadata,0,[]);

CD8 = extractChannels(36,metadata,0,[]);
hoechst = extractChannels(1,metadata,0,[]);
CD11c = extractChannels(43,metadata,0,[]);
CD206 = extractChannels(22,metadata,0,[]);
CD4 = extractChannels(26,metadata,0,[]);
CD103 = extractChannels(49,metadata,0,[]);

hoechstrs = imresize3(hoechst,0.125);
Ilog =filterLoGND(hoechstrs,2);
Imax = imregionalmin(Ilog);
thresh =  thresholdOtsu(Ilog(Imax));
nuclei  = (Ilog<thresh).*Imax;
nucleiEnlarged = imdilate(nuclei,strel('disk',3));

collagenrs = imresize3(collagen,0.125);
collagenth=imtophat(collagenrs,strel('disk',3));

threshold = thresholdOtsu(collagenth);
collagenmask = collagenrs>threshold;
Idist = bwdist(collagenmask);

CD4rs = imresize3(CD4,size(hoechstrs));
stats = regionprops3(bwlabeln(nucleiEnlarged,4),CD4rs,'MeanIntensity');
CD4threshold = thresholdOtsu(cat(1, stats.MeanIntensity));
idxCD4 = find([stats.MeanIntensity] > CD4threshold);
CD4nuclei = ismember(bwlabeln(nucleiEnlarged,4),idxCD4);  

CD103rs = imresize3(CD103,size(hoechstrs));
stats103 = regionprops3(bwlabeln(nucleiEnlarged),CD103rs,'MeanIntensity');
CD103threshold = multithresh(cat(1, stats103.MeanIntensity),2);
idxCD103 = find([stats103.MeanIntensity] > CD103threshold(1));
CD103nuclei = ismember(bwlabeln(nucleiEnlarged),idxCD103);  

CD8rs = imresize3(CD8,size(hoechstrs));
statsCD8 = regionprops3(bwlabeln(nucleiEnlarged),CD8rs,'MeanIntensity');
CD8threshold = multithresh(cat(1, statsCD8.MeanIntensity),2);
idxCD8 = find([statsCD8.MeanIntensity] > CD8threshold(1) & [stats103.MeanIntensity] < CD103threshold(1));
CD8nuclei = ismember(bwlabeln(nucleiEnlarged),idxCD8);  


CD11crs = imresize3(CD11c,size(hoechstrs));
stats = regionprops3(bwlabeln(nucleiEnlarged),CD11crs,'MeanIntensity');
CD11cthreshold = thresholdOtsu(cat(1, stats.MeanIntensity));
idxCD11c = find([stats.MeanIntensity] > CD11cthreshold);
CD11cnuclei = ismember(bwlabeln(nucleiEnlarged),idxCD11c);  

CD206rs = imresize3(CD206,size(hoechstrs));
stats = regionprops3(bwlabeln(nucleiEnlarged),CD206rs,'MeanIntensity');
CD206threshold = thresholdOtsu(cat(1, stats.MeanIntensity));
idxCD206 = find([stats.MeanIntensity] > CD206threshold);
CD206nuclei = ismember(bwlabeln(nucleiEnlarged),idxCD206);  


CD4dist = sum(sum(sum(Idist(CD4nuclei))))/sum(sum(sum(CD4nuclei)));
CD4diststd = std(Idist(CD4nuclei))/sqrt(numel(idxCD4));
CD8dist = sum(sum(sum(Idist(CD8nuclei))))/sum(sum(sum(CD8nuclei)));
CD8diststd = std(Idist(CD8nuclei))/sqrt(numel(idxCD8));
CD103dist = sum(sum(sum(Idist(CD103nuclei))))/sum(sum(sum(CD103nuclei)));
CD103diststd = std(Idist(CD103nuclei))/sqrt(numel(idxCD103));
CD11cdist = sum(sum(sum(Idist(CD11cnuclei))))/sum(sum(sum(CD11cnuclei)));
CD11cdiststd = std(Idist(CD11cnuclei))/sqrt(numel(idxCD11c));
CD206dist = sum(sum(sum(Idist(CD206nuclei))))/sum(sum(sum(CD206nuclei)));
CD206diststd = std(Idist(CD206nuclei))/sqrt(numel(idxCD206));

X = categorical({'CD4', 'CD8', 'CD103', 'CD206', 'CD11c'});
Y=[CD4dist,CD8dist,CD103dist, CD206dist,CD11cdist];
err = [CD4diststd,CD8diststd,CD103diststd, CD206diststd,CD11cdiststd];
errlow = Y-err;
errhigh = Y+err;
h=bar(X,Y,'FaceColor','flat')
h.CData(1,:) = [1 0 0];
h.CData(2,:) = [1 0 1];
h.CData(3,:) = [1 1 0];
h.CData(4,:) = [0 1 0];
h.CData(5,:) = [0 0 1];

hold on
h=errorbar(X,Y,err)
h.Color = [0 0 0];                            
h.LineStyle = 'none';  
%CD4 and CD206 are closer to collagen fibers



X = categorical({'CD103', 'CD8'});
Y=[CD103dist,CD8dist];
err = [CD103diststd,CD8diststd];

