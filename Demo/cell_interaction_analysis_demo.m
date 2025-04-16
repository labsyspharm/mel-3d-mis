channelA = extractChannels(1,bfGetReader('Figure 6.tif'),0,[]); %CD4
channelB = extractChannels(2,bfGetReader('Figure 6.tif'),0,[]); %CD8

%% tight interaction (between cell 2 and 3)
pixelSize = 0.14;
planeZ = 1;
xstart = 312;
xend = 338;
ystart = 349;
yend = 363;
rot_angle=0;
%% loose interaction (between cell 4 and 5)
pixelSize = 0.14;
planeZ = 4;
xstart = 425;
xend = 445;
ystart = 456;
yend = 458;
rot_angle=35;

%% analysis starts here
length = xend-xstart+1;
channelArotate = imrotate(channelA,rot_angle);
channelBrotate = imrotate(channelB,rot_angle);

test=[];
for i=ystart:yend
    test=cat(1,test,channelArotate(i,xstart:xend,planeZ));
end
test=mean(test,1);

test1=[];
for i=ystart:yend
    test1=cat(1,test1,channelBrotate(i,xstart:xend,planeZ));
end
test1=mean(test1,1);

newcolors = [0 0.7 0.7 ; 0.7 0 0.7 ; 0 0.7 0.7;  1 0 0; 0 0 0 ; 0.7 0 0.7;  1 0 0; 0 0 0  ];
colororder(newcolors)
plot((1:length)*pixelSize,test/max(test(:)),(1:length)*pixelSize,test1/max(test1(:)))

[f,gof]=fit((1:length)',test'/max(test(:)),'poly9');
hold on   

xf = 1:length;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*pixelSize,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:length);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
maxF=0;
for i=1:numel(unique_roots)
    if f(unique_roots(i))>maxF
        maxF = f(unique_roots(i));
        idx = i;
    end
end
plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')

[f,gof]=fit((1:length)',test1'/max(test1(:)),'poly9');
hold on   

xf = 1:length;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*pixelSize,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:length);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
maxF=0;
for i=1:numel(unique_roots)
    if f(unique_roots(i))>maxF
        maxF = f(unique_roots(i));
        idx = i;
    end
end
plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
legend ('cell1','cell2')

