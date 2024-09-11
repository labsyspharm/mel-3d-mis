metadata =bfGetReader(['V:\cycif-techdev\ClarenceLSM980\20microntonsil\tonsila 16bit bgsub.tif']);
allI = [];
for iChan=1:metadata.getSizeC
I = extractChannels(iChan,metadata,0,[]);
maxI = max(I,[],3);
allI= cat(5,allI,maxI);
disp(['Finished channel ' num2str(iChan)])
end
tiffwriteimj(allI,'D:\tonsila_MIP.tif')