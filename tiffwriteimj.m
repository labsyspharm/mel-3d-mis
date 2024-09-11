function tiffwriteimj(varargin)
% tiffwriteimj  Save a 5D matrix into ImageJ compatible TIFF
%
%   tiffwriteimj(I,outputPath,dimensionOrder,varargin)
%
%    tiffwriteimj(I, outputPath) writes the input 5D matrix into a new file
%    specified by outputPath.
%
%    tiffwriteimj(I, outputPath, dimensionOrder) specifies the dimension order of
%    the input matrix. Default valuse is YXZCT... (matlab puts Y first)
%
%    tiffwriteimj(I, outputPath, 'Compression', compression) specifies the
%    compression to use when writing the OME-TIFF file.
%
%    tiffwriteimj(I, outputPath, 'BigTiff', true) allows to save the file using
%    64-bit offsets
%
%    tiffwriteimj(I, outputPath, 'metadata', metadata) allows to pull
%    metadata from an OME-XML metadata object or a metadata struct:
%    A metadata struct can have the following fields:
%    'unit' (units of pixel dimensions provided)
%    'PixelSizeX'
%    'PixelSizeY'
%    'Zstep'
%    ... more to be added later ...
%
%    Examples:
%
%        tiffwriteimj(zeros(100, 100), outputPath)
%        tiffwriteimj(zeros(100, 100, 2, 3, 4), outputPath)
%        tiffwriteimj(zeros(100, 100, 20), outputPath, 'dimensionOrder', 'XYTZC')
%        tiffwriteimj(zeros(100, 100), outputPath, 'Compression', 'LZW')
%        tiffwriteimj(zeros(100, 100), outputPath, 'BigTiff', true)
%        tiffwriteimj(zeros(100, 100), outputPath, 'metadata', metadata)
%
% Talley Lambert
% 4/2016

    % Input check
    ip = inputParser;
    ip.addRequired('I', @isnumeric);
    ip.addRequired('outputPath', @ischar);
    % default dimension order for matlab is YXZ ...
    ip.addOptional('dimensionOrder', 'YXZCT', @(x) ismember(x, {'XYZCT','XYZTC','XYCZT','XYTZC','XYCTZ','XYTCZ','YXZCT','YXZTC','YXCZT','YXTZC','YXCTZ','YXTCZ'}));
    ip.addParameter('metadata', [], @(x) isa(x, 'loci.formats.ome.OMEXMLMetadata') || isstruct(x));
    ip.addParameter('Compression', '',  @(x) ismember(x, fieldnames(Tiff.Compression)));
    ip.addParameter('BigTiff', false , @islogical);
    ip.parse(varargin{:});
    
    % rearrange the stack to be in order XYCZT (imageJ default)
    % this was originally XYZCT...
    % and get the dimension sizes
    permuteOrder=zeros(1,5);
    for q=1:5
        switch ip.Results.dimensionOrder(q)
            case {'X','x'}
                permuteOrder(1)=q;
                x = size(ip.Results.I,q);
            case {'Y','y'}
                permuteOrder(2)=q;
                y = size(ip.Results.I,q);
            case {'C','c'}
                permuteOrder(3)=q;
                c = size(ip.Results.I,q);
            case {'Z','z'}
                permuteOrder(4)=q;
                z = size(ip.Results.I,q);
            case {'T','t'}
                permuteOrder(5)=q;
                t = size(ip.Results.I,q);
        end
    end
    array=permute(ip.Results.I,permuteOrder);

    % begin building Tiff Tags
    tagstruct.ImageLength = y;
    tagstruct.ImageWidth = x;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.ResolutionUnit=Tiff.ResolutionUnit.Centimeter;
    %tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.Software = 'MATLAB';
    tagstruct.XResolution = 1;
    tagstruct.YResolution = 1;

    
    % Create metadata struct based on array data
    md.ImageJ='1.50i';
    md.images=c*z*t;
    md.channels=c;
    md.slices=z;
    md.frames=t;
    md.spacing=1;
    md.finterval=0; % this is unused at the moment... could store frame interval?
    md.hyperstack='true';
    md.mode='grayscale';
    md.loop='false';
    md.unit='micron';

    if isa(ip.Results.metadata, 'loci.formats.ome.OMEXMLMetadata')
        
        scale=double(ip.Results.metadata.getPixelsPhysicalSizeX(0).unit.getScaleFactor);
        tagstruct.XResolution = (1e-06/scale)/double(ip.Results.metadata.getPixelsPhysicalSizeX(0).value);
        tagstruct.YResolution = (1e-06/scale)/double(ip.Results.metadata.getPixelsPhysicalSizeY(0).value);
        
        md.spacing=double(ip.Results.metadata.getPixelsPhysicalSizeZ(0).value);
        
        if ip.Results.metadata.getPixelsSizeX(0)~=x;
            warning('OME Metadata provided does not agree with X dimension of input array')
        end
        if ip.Results.metadata.getPixelsSizeY(0)~=y;
            warning('OME Metadata provided does not agree with Y dimension of input array')
        end
        if ip.Results.metadata.getPixelsSizeZ(0)~=z;
            warning('OME Metadata provided does not agree with Z dimension of input array')
        end
        if ip.Results.metadata.getPixelsSizeC(0)~=c;
            warning('OME Metadata provided does not agree with C dimension of input array')
        end
        if ip.Results.metadata.getPixelsSizeT(0)~=t;
            warning('OME Metadata provided does not agree with T dimension of input array')
        end

        % could do lots more with OMEXML data here....
        
    elseif isstruct(ip.Results.metadata)

        % assume microns unless unit is specified
        % but then convert to microns
        if isfield(ip.Results.metadata,'unit')
            switch ip.Results.metadata.unit
                case {'micrometer','micron','um','µm'}
                    scale=1e-6;
                case {'centimeter','cm'}
                    scale=1e-2;
                case {'millimeter','mm'}
                    scale=1e-3;
                case {'nanometer','nm'}
                    scale=1e-9;
                case {'meter','m'}
                    scale=1;
            end
        else
            scale=1e-06;
        end
            
        tagstruct.XResolution = (1e-06/scale)/ip.Results.metadata.PixelSizeX;
        tagstruct.YResolution = (1e-06/scale)/ip.Results.metadata.PixelSizeY;
        md.spacing=ip.Results.metadata.Zstep/(1e-06/scale);
    end

    
    % type detection
    switch class(array)
        case {'uint8'}
            tagstruct.BitsPerSample = 8;
            tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        case {'uint16'}
            tagstruct.BitsPerSample = 16;
            tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        case {'uint32'}
            tagstruct.BitsPerSample = 32;
            tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        case {'int8'}
            array=uint8(array);
            tagstruct.BitsPerSample = 8;
            tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        case {'int16'}
            array=uint16(array);
            tagstruct.BitsPerSample = 16;
            tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        case {'single'}
            tagstruct.BitsPerSample = 32;
            tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        case 'double'
            array=single(array);
            tagstruct.BitsPerSample = 32;
            tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    end
    
    
    % Create ImageWriter
    if ip.Results.BigTiff
        outtiff=Tiff(ip.Results.outputPath, 'w8');
    else
        % check whether data is greater than 4GB and force big Tiff
        % probably should allow some user control over this
        if numel(array)*tagstruct.BitsPerSample/8 > 4000000000
            outtiff=Tiff(ip.Results.outputPath, 'w8');
        else
            outtiff=Tiff(ip.Results.outputPath, 'w');
        end
    end
    
    if ~isempty(ip.Results.Compression)
        tagstruct.Compression = Tiff.Compression.(ip.Results.Compression);
    end
    
    
    % It'd be nice to have a seperate min/max display range for each
    % channel but it doesn't appear to be possible in ImageJ with just the
    % 270 tiff tag alone... 
    if min(array(:)) < 0
        md.min=0;
    else
        md.min=min(array(:));
    end
    md.max=max(array(:));
    
    
    % turn metadata into string
    fields = fieldnames(md);
    celery = cell(1,numel(fields));
    for f = 1:numel(fields)
        celery{f}=(strcat(fields{f},'=',num2str(md.(fields{f}))));
    end
    metadata=strjoin(celery,char(13)); % char(13) is ASCII newline...
    tagstruct.ImageDescription = metadata;
    
    % write stack
    zctCoord = [c z t];
    for index = 1 : (c*z*t)
        [i, j, k] = ind2sub(zctCoord, index);
        outtiff.setTag(tagstruct);
        plane = array(:, :, i, j, k)';
        outtiff.write(plane);
        outtiff.writeDirectory();
    end
    outtiff.close()
    
end