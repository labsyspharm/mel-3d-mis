metadata =bfGetReader(['D:\F8iia_MIS_MIP.tif']);
channels = [{58 14 4 10 3 0};
    {15 7 0 2 20 0};
    {23 26 36 11 0 16};
    {22 19 30 31 35 0};
    {34 42 39 38 32 0};
    {44 51 49 46 47 0};
    {54 63 57 56 52 0};
    {4 66 61 27 69 40};
    {59 68 8 28 48 0};
    {0 0 43 55 64 70}];

allI=[];
for iPanel=1:10
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = I.*0;
        elseif iChan==70
            I = imadjust(extractChannels(iChan,metadata,0,[]));
            I=imresize(I,0.25,'bicubic');
            I=imrotate(I,-25);
            I=I(700:1700,:);
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0 0.9]);
             I=imresize(I,0.25,'bicubic');
            I=imrotate(I,-25);
            I=I(700:1700,:);
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\F8iia_MIS_MIP_gallery.tif')



metadata =bfGetReader(['D:\MIP_Dataset1-LSP13626-invasive_margin.tif']);
channels = [{58 14 4 10 3 0};
    {15 7 0 2 20 0};
    {23 26 36 11 0 16};
    {22 19 30 31 35 0};
    {34 42 39 38 32 0};
    {44 51 49 46 47 0};
    {54 63 57 56 52 0};
    {4 66 61 27 69 40};
    {59 68 8 28 48 0};
    {0 0 43 55 64 0}];
allI=[];
for iPanel=1:10
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = I.*0;
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0 0.9]);
             I=imresize(I,0.25,'bicubic');
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\F8iic_IM_MIP_gallery.tif')



metadata =bfGetReader(['D:\LSP17378_GBM_MIP.tif']);
channels = [{2 36 28 6 1};
           {7 8 10 11 0};
           {12 14 15 16 0};
           {0 19 20 22 0};
           {23 40 26 27 0};
           {30 31 32 34 3}];

allI=[];
for iPanel=1:6
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = I.*0;
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0 0.9]);
             I=imresize(I,0.25,'bicubic');
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\GBM_MIP_gallery.tif')




metadata =bfGetReader(['U:\cycif-techdev\confocal\40xMET15microns\outputb\MAX_cycle1_3D_MERGED.tif']);
channels = [{1 2 3};
    {4 5 6};
    {7 8 9};
    {10 11 12};
    {13 14 15};
    {16 17 18};
];

allI=[];
for iPanel=1:6
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = I.*0;
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0 0.99]);
%              I=imresize(I,0.5,'bicubic');
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\lungMET_MIS_MIP_gallery.tif')



metadata =bfGetReader(['V:\cycif-techdev\ClarenceLSM980\20micronF8\MIP_F8a 16bit bgsub.tif']);
channels = [{0 2 3 4 1};
    {52 8 9 11 0};
    {13 14 16 17 0} ;
    {19 21 22 23 0};
    {26 27 28 29 0};
    {31 32 33 34 0};
    {37 38 39 41 0};
    {43 44 46 47 0};
    {6 12 24 36 0};
    {42 48 50 51 0};
];

allI=[];
for iPanel=1:10
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = I.*0;
        elseif iChan ==1
            I = extractChannels(iChan,metadata,0,[]);
             I=imresize(I,0.25,'bicubic');
            I=imrotate(I,-25);
            I=I(750:1700,:);
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0 0.9]);
             I=imresize(I,0.25,'bicubic');
            I=imrotate(I,-25);
            I=I(750:1700,:);
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\F8a_MIS_MIP_gallery.tif')


metadata =bfGetReader(['V:\cycif-techdev\ClarenceLSM980\20micronF8\MIP_F8c 16bit bgsub.tif']);
channels = [{1 2 3 4 0};
    {0 52 8 9 11 };
    {0 13 14 16 17 } ;
    {0 19 21 22 23 };
    {0 26 27 28 29 };
    {0 31 32 33 34 };
    {0 37 38 39 41 };
    {0 43 44 46 47 };
    {0 6 12 24 36 };
    {0 42 48 50 51 };
];

allI=[];
for iPanel=1:10
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = I.*0;
              
        elseif iChan ==1
            I = extractChannels(iChan,metadata,0,[]);
             I=imresize(I,0.25,'bicubic');
           I=imrotate(I,-90);
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0 0.9]);
             I=imresize(I,0.25,'bicubic');
           I=imrotate(I,-90);
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\F8c_MIS_MIP_gallery.tif')



metadata =bfGetReader(['V:\cycif-techdev\ClarenceLSM980\20microntonsil\MIP_tonsilc 16bit bgsub.tif']);
channels = [{20 4 5 6};
    {0 16 9 10};
    {29 11 13 14} ;
    {0 1 18 19};
    {0 22 23 24};
    {0 27 28 31};
    {0 32 33 34};
    {7 8 15 26};
];

allI=[];
for iPanel=1:8
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = I.*0;
            
          
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0 0.9]);
             I=imresize(I,0.25,'bicubic');
           
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\tonsilc_MIP_gallery.tif')



metadata =bfGetReader(['U:\cycif-production\110-BRCA-Mutant-Ovarian-Precursors\John_STIC_3D\ovarianTR4_MIP.tif']);
channels = [{0 6 3 4 1};
    {7 8 10 11 0};
    {12 14 15 16 0} ;
    {19 20 22 23 0};
    {27 28 30 31 0 };
    {0 48 36 38 0 };
    {40 42 43 44 0 };
    {47 50 51 52 0 };
    {11 18 24 32 39 };
    {46 2 0 0 0 };
];

allI=[];
for iPanel=1:10
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = I.*0;
            
          
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0.2 0.99]);
             I=imresize(I,0.25,'bicubic');
           
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\STICtr4_MIP_gallery.tif')



metadata =bfGetReader(['U:\cycif-production\110-BRCA-Mutant-Ovarian-Precursors\John_STIC_3D\ovarianTR5_MIP.tif']);
channels = [{0 6 3 4 1};
    {7 8 10 11 0};
    {12 14 15 16 0} ;
    {19 20 22 23 0};
    {27 28 30 31 0 };
    {0 48 36 38 0 };
    {40 42 43 44 0 };
    {47 50 51 52 0 };
    {11 18 24 32 39 };
    {46 2 0 0 0 };
];

allI=[];
for iPanel=1:10
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = imadjust(extractChannels(1,metadata,0,[]),[],[0.2 0.99]);
            I = I.*0;
          I=imresize(I,0.25,'bicubic');
          I=imrotate(I,-90);
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0.2 0.99]);
             I=imresize(I,0.25,'bicubic');
           I=imrotate(I,-90);
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\STICtr5_MIP_gallery.tif')

metadata =bfGetReader(['U:\cycif-production\110-BRCA-Mutant-Ovarian-Precursors\John_STIC_3D\ovarianTR3_MIP.tif']);
channels = [{0 6 3 4 1};
    {7 8 10 11 0};
    {12 14 15 16 0} ;
    {19 20 22 23 0};
    {27 28 30 31 0 };
    {0 48 36 38 0 };
    {40 42 43 44 0 };
    {47 50 51 52 0 };
    {11 18 24 32 39 };
    {46 2 0 0 0 };
];

allI=[];
for iPanel=1:10
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = imadjust(extractChannels(1,metadata,0,[]),[],[0.2 0.99]);
            I = I.*0;
          I=imresize(I,0.25,'bicubic');
          I=imrotate(I,-90);
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0.2 0.99]);
             I=imresize(I,0.25,'bicubic');
           I=imrotate(I,-90);
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\STICtr3_MIP_gallery.tif')


metadata =bfGetReader(['U:\cycif-techdev\Clarence_LSM980\20umCHTN\MIP_CHTN 16bit bgsub.tif']);
channels = [{1 2 3 4 5};
    {7 9 11 12 14};
    {15 16 17 20 21 } ;
    {23 24 25 26 0};
    {28 29 30 31 33 };
    {34 35 37 38 39};
    {41 42 43 45 46 };
    {47 48 50 51 52};
    {55 56 57 59 60 };
    {61 64 65 0 0 };
];

allI=[];
for iPanel=1:10
    panel = [];
    for iChan=cat(2,channels{iPanel,:})
        if iChan ==0
            I = imadjust(extractChannels(1,metadata,0,[]),[],[0.2 0.99]);
            I = I.*0;
          I=imresize(I,0.25,'bicubic');
        elseif iChan ==1
          I = extractChannels(1,metadata,0,[]);
          I=imresize(I,0.25,'bicubic');
        else
            I = imadjust(extractChannels(iChan,metadata,0,[]),[],[0.2 0.99]);
             I=imresize(I,0.25,'bicubic');
        end
        panel= cat(4,panel,I);
    end
        disp(['Finished channel ' num2str(iChan)])
allI = cat(5,allI,panel);
end
tiffwriteimj(allI,'D:\CHTN_MIP_gallery.tif')


