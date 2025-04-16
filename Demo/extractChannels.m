

function V= extractChannels(k,metadata,angle,crop)
nZ = metadata.getSizeZ;
nr = metadata.getSizeY;
nc = metadata.getSizeX;
dType = metadata.getBitsPerPixel;
nChan = metadata.getSizeC;
totalPlanes = nZ*nChan;

if ~isempty(crop)
    nr = crop(4);
    nc = crop(3);
end

V = zeros(nr,nc,1,'uint16');
    for iZ = 1:nZ
        disp(['Reading ' int2str(k) ' of ' int2str(totalPlanes)])
        if dType == 8
            if isempty(crop)
                V(:,:,iZ) = imrotate(bfGetPlane(metadata,k),angle);
            else
                V(:,:,iZ) = imrotate(bfGetPlane(metadata,k,crop(1),crop(2),crop(3),crop(4)),angle);
            end
        elseif dType == 16
%             V(:,:,iZ) = imcrop(imrotate(uint8(double(bfGetPlane(metadata,k))/65535*255),angle),[xStart,yStart,width,height]);
%             V(:,:,iZ) = uint8(double(bfGetPlane(metadata,k))/65535*255);
            if isempty(crop)
                V(:,:,iZ) = imrotate(bfGetPlane(metadata,k),angle);
%                             V(:,:,iZ) = uint8(double(bfGetPlane(metadata,k))/65535*255);
            else
%                 V(:,:,iZ) = imrotate(uint8(double(bfGetPlane(metadata,k,crop(1),crop(2),crop(3),crop(4)))/65535*255),angle);
                V(:,:,iZ) = imrotate(bfGetPlane(metadata,k,crop(1),crop(2),crop(3),crop(4)),angle);
            end
        else
            print('Unexpected data type detected')
            return
        end
        k=k+nChan;
    end
end