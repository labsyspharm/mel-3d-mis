% function mashlar(fileNum)
bfCheckJavaPath()
% 
% 
iROI=1;
fixedCycle=1;
DAPIchan=1;
mainPath = 'U:\cycif-techdev\Clarence_LSM980\20umF8-IIprocessed';
outputPath = ['F:\'];
mkdir(outputPath)
ROI = {'a';'b';'c';'d';};
realfixedposX = 2700;
realfixedposY = 0;
offsetX = 5.3264e+05; %MIS
offsetY = 2.6509e+05;
% offsetX = 4.9953e+05; %IM
% offsetY = 3.5003e+05;
scalingFactor =0.25;
scalingFactorZ =0.25;
driftCorrect=1;
% 
% 
cycleFolder = dir([mainPath filesep 'F8' ROI{iROI} '-cycle*Processing.czi']);
cycleFolder = natsortfiles(cycleFolder);
disp(['Processing ROI ' ROI{iROI}])

% find largest Z stack
largestZ=169;
% for iCycle = 1:numel(cycleFolder)
%     disp(['Inspecting ' cycleFolder(iCycle).folder filesep cycleFolder(iCycle).name])
%     metadata = loci.formats.Memoizer(bfGetReader());
%     metadata.setId([cycleFolder(iCycle).folder filesep cycleFolder(iCycle).name])
%     if metadata.getSizeZ > largestZ
%         largestZ =  metadata.getSizeZ;
%     end
% end

%load fixed cycle
r = loci.formats.Memoizer(bfGetReader());
options = javaObject('loci.formats.in.DynamicMetadataOptions');
options.setBoolean(java.lang.String('zeissczi.autostitch'), java.lang.Boolean('TRUE'));
r.setMetadataOptions(options);
r.setId('U:\cycif-techdev\Clarence_LSM980\20umF8-IIprocessed\F8a-cycle1-LSM Plus Processing-Stitching.czi');
omeMetadata = r.getMetadataStore();
disp(num2str(r.getSeriesCount()))
fixedposX = omeMetadata.getPlanePositionX(1,1).value.doubleValue/omeMetadata.getPixelsPhysicalSizeX(0).value.doubleValue;
fixedposY = omeMetadata.getPlanePositionY(1,1).value.doubleValue/omeMetadata.getPixelsPhysicalSizeY(0).value.doubleValue;

volReshaped =  volumeRead8bitConvert('U:\cycif-techdev\Clarence_LSM980\20umF8-IIprocessed\F8a-cycle1-LSM Plus Processing-Stitching.czi',0,0,[1 1 9585 4000])-0;
nChan = size(volReshaped,4);
disp('First cycle loaded')

%             volReshaped =  extractChannels(1,r,0);


% find best plane and insert cycle 1 into final array
sumstd =0;
sharpZ=0;
test = stdfilt(imresize(volReshaped(:,:,:,DAPIchan),0.25));
for iZ = 1:size(volReshaped(:,:,:,DAPIchan),3)
   sumstdplane =  sum(sum(test(:,:,iZ)));
   if sumstd <sumstdplane
       sumstd = sumstdplane;
       sharpZ = iZ;
   end
end
% volReg = cat(3,zeros(size(volReshaped,1),size(volReshaped,2),round(largestZ*1.15/2)-sharpZ, size(volReshaped,4),'uint8'), volReshaped);
% volReg = cat(3,volReg,zeros(size(volReshaped,1),size(volReshaped,2),round(largestZ*1.15/2)-(size(volReshaped,3)-sharpZ), size(volReshaped,4), 'uint8'));
volReg = cat(3,zeros(size(volReshaped,1),size(volReshaped,2),round(largestZ*1.15/2)-sharpZ, size(volReshaped,4),'uint16'), volReshaped);
volReg = cat(3,volReg,zeros(size(volReshaped,1),size(volReshaped,2),round(largestZ*1.15/2)-(size(volReshaped,3)-sharpZ), size(volReshaped,4), 'uint16'));
clear volReshaped

if fileNum==1
    metadata = createMinimalOMEXMLMetadata(volReg);
    pixelSize = ome.units.quantity.Length(java.lang.Double(0.140), ome.units.UNITS.MICROMETER);
    metadata.setPixelsPhysicalSizeX(pixelSize, 0);
    metadata.setPixelsPhysicalSizeY(pixelSize, 0);
    metadata.setPixelsPhysicalSizeZ(pixelSize, 0);
    pixelSizeZ = ome.units.quantity.Length(java.lang.Double(0.280), ome.units.UNITS.MICROMETER);
    metadata.setPixelsPhysicalSizeZ(pixelSizeZ, 0);
    tiffwriteimj(volReg,[outputPath filesep 'F8II' ROI{iROI} '_cycle1_MERGEDNoDemons.tif'],'metadata',metadata,'dimensionOrder', 'XYZCT')
    return
else
    disp(['Proceeding to cycle' int2str(fileNum)])
end

volReg = volReg(:,:,:,DAPIchan);
maxDAPIref = max(volReg,[],3);


%% register subsequent cycles
for iCycle=fileNum:fileNum%2:numel(cycleFolder)    
    metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
    options.setBoolean(java.lang.String('zeissczi.autostitch'), java.lang.Boolean('FALSE'));
    metadata.setMetadataOptions(options);
    metadata.setId([mainPath filesep cycleFolder(iCycle).name]);
    numChan = metadata.getSizeC;
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
%     final = zeros(size(volReg,1),size(volReg,2),size(volReg,3),numChan,'uint8');
    final = zeros(size(volReg,1),size(volReg,2),size(volReg,3),numChan,'uint16');
    if fileNum ==16 ||fileNum ==17
        bgoffset = 10050
    else
        bgoffset = 0
    end
        for iPlane = 1:metadata.getSeriesCount()
            tic
            [I,posX,posY] = bftileopen([mainPath filesep cycleFolder(iCycle).name],iPlane);
            posX = posX - offsetX;
            posY = posY - offsetY;
    
            vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
               for iZ = 1:sizeZ
                   for iC =  1:numChan
                        vol(:,:,iZ,iC) = I{iPlane,1}{(iZ-1)*numChan+iC} - bgoffset;
                   end
               end
            clear I
            DAPI = vol(:,:,:,1);
            maxI = max(DAPI,[],3);
            tileSize = 1330;
            ystart = posY - tileSize*0.3;
            xstart = posX - tileSize*0.3; 
            yend = posY + tileSize + tileSize*0.3;
            xend = posX + tileSize + tileSize*0.3;
    
            if ystart<1
                ystart=1;
            end
            if xstart<1
                xstart=1;
            end
            if yend > size(maxDAPIref,1)
                yend = size(maxDAPIref,1);
            end
            if xend > size(maxDAPIref,2)
                xend = size(maxDAPIref,2);
            end
    
            tile = maxDAPIref(ystart:yend,xstart:xend);
            tile3D = volReg(ystart:yend,xstart:xend,:);
    %         figure,imshow(tile,[])
    
          
            % load fixed target
            [optimizer, metric] = imregconfig('monomodal');
            if sum(sum(tile))>10000
%                 tform2D= imregcorr(imresize(uint8(double(maxI)/65535*255)-40,scalingFactor),imresize(tile,scalingFactor),'transformType','translation');
                    tform2D= imregcorr(imresize(maxI,scalingFactor),imresize(tile,scalingFactor),'transformType','translation');
                init2Dtransform = eye(4,4);       
                init2Dtransform(4,1)=tform2D.T(3,1)/scalingFactor;
                init2Dtransform(4,2)=tform2D.T(3,2)/scalingFactor;
%                 prereg = zeros(size(tile3D),'uint8');
                    prereg = zeros(size(tile3D),'uint16');
                    for iChan = 1:numChan
                         test=imwarp(vol(:,:,:,iChan),affine3d(init2Dtransform),'OutputView',imref3d(size(tile3D)),'interp','nearest');
%                          prereg(:,:,:,iChan) = uint8(double(test)/65535*255);
                            prereg(:,:,:,iChan) = uint16(test);
                    end
                clear test
    
                figure,imshowpair(max(prereg(:,:,:,1),[],3),max(tile,[],3))
                [optimizer, metric] = imregconfig('monomodal');
                tform = imregtform(imresize3(prereg(:,:,:,1),'Scale',[scalingFactor scalingFactor scalingFactorZ]),...
                    imresize3(tile3D,'Scale',[scalingFactor scalingFactor scalingFactorZ]),'translation', optimizer,metric);
    
                tform.T(4,1:2)=1;
                tform.T(4,3) =tform.T(4,3)/scalingFactorZ; 
                disp(num2str(tform.T(4,3)))
               
              

                for iChan = 1:numChan
                registered=imwarp(prereg(:,:,:,iChan),tform,'OutputView',imref3d(size(tile3D)),'interp','nearest');
%                 imshowpair(max(permute(registered,[2 3 1])-40,[],3),max(permute(tile3D,[2 3 1])-40,[],3))
                if iChan ==1 && driftCorrect==1
                    D = imregdemons(imresize3(registered-40,0.125),imhistmatchn(imresize3(tile3D,0.125),registered),'AccumulatedFieldSmoothing',1);
                    D=gather(D);
                end
                
                if driftCorrect==1
                    upsampledD = zeros([size(tile3D) 3],'int8');
                    for iZ = 1:size(D,4)
                        upsampledD(:,:,:,iZ) = 8*imresize3(D(:,:,:,iZ),size(tile3D),'nearest');
                    end
%                     registered = uint8(imwarp(registered,upsampledD,'interp','nearest'));
                registered = imwarp(registered,upsampledD,'interp','nearest');
                end
                background = final(ystart:yend,xstart:xend,:,iChan);
                final(ystart:yend,xstart:xend,:,iChan) = background + registered.*uint16(~(background>0));
                
                end
                clear registered prereg 


            else
                continue
            end
            toc
        end
        
    
    metadata = createMinimalOMEXMLMetadata(final);
    pixelSize = ome.units.quantity.Length(java.lang.Double(0.140), ome.units.UNITS.MICROMETER);
    metadata.setPixelsPhysicalSizeX(pixelSize, 0);
    metadata.setPixelsPhysicalSizeY(pixelSize, 0);
    metadata.setPixelsPhysicalSizeZ(pixelSize, 0);
    pixelSizeZ = ome.units.quantity.Length(java.lang.Double(0.280), ome.units.UNITS.MICROMETER);
    metadata.setPixelsPhysicalSizeZ(pixelSizeZ, 0);
    tiffwriteimj(final,[outputPath filesep 'F8II' ROI{iROI} '_cycle' int2str(iCycle) '_MERGEDNoDemons.tif'],'metadata',metadata,'dimensionOrder', 'XYZCT')
    
end
disp(['Completed ' ROI{iROI} ])

% end
