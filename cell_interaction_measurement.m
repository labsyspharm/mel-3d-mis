%%%%%%%%%%%%%% - IM Direct interactions PD1 PDL1 on immune cells
PDL1 = extractChannels(2,bfGetReader('D:\IM_direct interaction PDL1CD4PD1CD11c.tif'),0,[]);
PD1 = extractChannels(4,bfGetReader('D:\IM_direct interaction PDL1CD4PD1CD11c.tif'),0,[]);
PDL1rotate = imrotate(PDL1,120);
PD1rotate = imrotate(PD1,120);
test=[];
for i=276:284
    test=cat(1,test,PDL1rotate(i,333:353,49));
end
test=mean(test,1);

test1=[];
for i=276:284
    test1=cat(1,test1,PD1rotate(i,333:353,49));
end
test1=mean(test1,1);
newcolors = [0 1 0 ; 1 0 1 ; 0 1 0 ;  1 0 0; 0 0 0; 1 0 1 ; 1 0 0; 0 0 0];
colororder(newcolors)

[f,gof]=fit((1:21)',test(1:21)'/max(test(:)),'poly9');
colororder(newcolors)
plot((1:21)*0.14, test(1:21)/max(test(:)), (1:21)*0.14, test1(1:21)/max(test1(:))) 
hold on   
xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')



[f,gof]=fit((1:21)',test1(1:21)'/max(test1(:)),'poly9');
colororder(newcolors)
hold on   
xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')

legend off

%%%%%%%%%%%%%% - IM Direct interactions PD1 PDL1 on tumor
PDL1 = extractChannels(1,bfGetReader('D:\MIS_PD1PDL1-tumor-tightinteraction.tif'),0,[]);
PD1 = extractChannels(2,bfGetReader('D:\MIS_PD1PDL1-tumor-tightinteraction.tif'),0,[]);
PDL1rotate = imrotate(PDL1,40);
PD1rotate = imrotate(PD1,40);
test=[];
for i=301:304
    test=cat(1,test,PDL1rotate(i,310:330,86));
end
test=mean(test,1);

test1=[];
for i=301:304
    test1=cat(1,test1,PD1rotate(i,310:330,86));
end
test1=mean(test1,1);
newcolors = [0 1 0 ; 1 0 1 ; 0 1 0 ;  1 0 0; 0 0 0; 1 0 1 ; 1 0 0; 0 0 0];
colororder(newcolors)

[f,gof]=fit((1:21)',test(1:21)'/max(test(:)),'poly9');
colororder(newcolors)
plot((1:21)*0.14, test(1:21)/max(test(:)), (1:21)*0.14, test1(1:21)/max(test1(:))) 
hold on   
xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')



[f,gof]=fit((1:21)',test1(1:21)'/max(test1(:)),'poly9');
hold on   
xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CD4 = extractChannels(1,bfGetReader('D:\MIS_CD4CD3CD11c_tightinteraction.tif'),0,[]);
CD3 = extractChannels(2,bfGetReader('D:\MIS_CD4CD3CD11c_tightinteraction.tif'),0,[]);
CD11c = extractChannels(3,bfGetReader('D:\MIS_CD4CD3CD11c_tightinteraction.tif'),0,[]);

%CD4-CD11c
test=[];
for i=243:261
    test=cat(1,test,CD4(i,183:205,66));
end
test=mean(test,1);

test1=[];
for i=243:261
    test1=cat(1,test1,CD11c(i,183:205,66));
end
test1=mean(test1,1);
newcolors = [1 0 1 ; 0.7 0.7 0];
colororder(newcolors)
plot((1:23)*0.14,[test/max(test(:)) ; test1/max(test1(:))])


%CD4CD3+CD11c type 2 interaction
CD4rotate = imrotate(CD4,-55);
CD3rotate = imrotate(CD3,-55);
CD11crotate = imrotate(CD11c,-55);
test=[];
for i=334:336
    test=cat(1,test,CD3rotate(i,334:358,63));
end
test=mean(test,1);

test1=[];
for i=334:336
    test1=cat(1,test1,CD11crotate(i,334:358,63));
end
test1=mean(test1,1);

test2=[];
for i=334:336
    test2=cat(1,test2,CD4rotate(i,334:358,63));
end
test2=mean(test2,1);

newcolors = [0 1 1 ; 1 0 1 ; 1 1 0; 0 1 1; 1 0 0; 0 0 0; 1 0 1; 1 0 0; 0 0 0; 1 1 0 ; 1 0 0; 0 0 0];
colororder(newcolors)
plot((1:25)*0.14,test/max(test(:)),(1:25)*0.14,test1/max(test1(:)),(1:25)*0.14,test2/max(test2(:)))

[f,gof]=fit((1:25)',test(1:25)'/max(test(:)),'poly9');
colororder(newcolors)
hold on   
xf = 1:25;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:25);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(3)*0.14,f(unique_roots(3)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(3) unique_roots(3)]*0.14, [0 1],'-')


[f,gof]=fit((1:25)',test1(1:25)'/max(test1(:)),'poly9');
hold on   
xf = 1:25;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:25);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(5)*0.14,f(unique_roots(5)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(5) unique_roots(5)]*0.14, [0 1],'-')


[f,gof]=fit((1:25)',test2(1:25)'/max(test2(:)),'poly9');
hold on   
xf = 1:25;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:25);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test2)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test2)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')


%%%%%%%%%%%%%%
CD4 = extractChannels(3,bfGetReader('D:\MIS_5 cells immune tightness.tif'),0,[]);
CD8 = extractChannels(6,bfGetReader('D:\MIS_5 cells immune tightness.tif'),0,[]);
MART1 = extractChannels(2,bfGetReader('D:\MIS_5 cells immune tightness.tif'),0,[]);

%% tighter junction
CD4rotate = imrotate(CD4,90);
MART1rotate = imrotate(MART1,90);

test=[];
for i=382:387
    test=cat(1,test,CD4rotate(i,321:341,79));
end
test=mean(test,1);

test1=[];
for i=382:387
    test1=cat(1,test1,MART1rotate(i,321:341,79));
end
test1=mean(test1,1);

newcolors = [0 0.7 0.7 ; 0 0.8 0 ; 0 0.7 0.7; 0 0 0; 0 0 0;0 0.8 0; 0 0 0];
colororder(newcolors)
plot((1:21)*0.14,test/max(test(:)),(1:21)*0.14,test1/max(test1(:)))


[f,gof]=fit((1:21)',test'/max(test(:)),'poly9');
hold on   

xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 
 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')

[f,gof]=fit((1:21)',test1'/max(test1(:)),'poly9');
hold on   

xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')
legend off


%% looser junction
CD4rotate = imrotate(CD4,45);
MART1rotate = imrotate(MART1,45);

test=[];
for i=530:535
    test=cat(1,test,CD4rotate(i,435:455,79));
end
test=mean(test,1);

test1=[];
for i=530:535
    test1=cat(1,test1,MART1rotate(i,435:455,79));
end
test1=mean(test1,1);

newcolors = [0 0.7 0.7 ; 0 0.8 0 ; 0 0.7 0.7; 0 0 0; 0 0 0;0 0.8 0; 0 0 0];
colororder(newcolors)
plot((1:21)*0.14,test/max(test(:)),(1:21)*0.14,test1/max(test1(:)))


[f,gof]=fit((1:21)',test'/max(test(:)),'poly9');
hold on   

xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 
 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(5)*0.14,f(unique_roots(5)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(5) unique_roots(5)]*0.14, [0 1],'-')

[f,gof]=fit((1:21)',test1'/max(test1(:)),'poly9');
hold on   

xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(1)*0.14,f(unique_roots(1)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(1) unique_roots(1)]*0.14, [0 1],'-')
legend off




% %% MART1 CD4 type 2 interaction
% MART1permute = permute(MART1,[3 2 1]);
% CD4permute = permute(CD4,[3 2 1]);
% test=[];
% for i=68:72
%     test=cat(1,test,CD4permute(i,257:287,352));
% end
% test=mean(test,1);
% 
% test1=[];
% for i=68:72
%     test1=cat(1,test1,MART1permute(i,257:287,352));
% end
% test1=mean(test1,1);
% 
% newcolors = [0 0.7 0.7 ; 0 0.8 0 ; 1 0 1];
% colororder(newcolors)
% plot((1:31)*0.14,test/max(test(:)),(1:31)*0.14,test1/max(test1(:)))
% 
% 
% lineintegral = max(cat(1,test/max(test(:)), test1/max(test1(:))));
% [f,gof]=fit((1:31)',lineintegral','poly9');
% newcolors = [0 0.7 0.7 ; 0 0.8 0 ; 0 0 0; 0 0 0; 0 0 0;0 0 0; 0 0 0];
% colororder(newcolors)
% plot((1:28)*0.14, test(1:28)/max(test(:)), (1:28)*0.14, test1(1:28)/max(test1(:))) 
% hold on   
% 
% xf = 1:28;
% yf=[];
% for i=xf
%     yf=cat(2,yf,f(xf(i)));
% end
% plot(xf*0.14,yf,'--')
% 
%  [fx fxx]= differentiate(f,1:0.5:31);
%               fdxHan=@(x)(differentiate(f,x));
%  dxroots =[];
%               for i = 1:numel(lineintegral)
%                 dxroots(i)=fzero(fdxHan, i);
%                 if dxroots(i) <0 || dxroots(i) > numel(lineintegral)
%                     dxroots(i)=NaN;
%                 end
%               end
% unique_roots=uniquetol(dxroots,0.001);
% plot(unique_roots(2)*0.14,f(unique_roots(2)),'x','MarkerEdgeColor','red','MarkerSize',12)
% plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
% plot([unique_roots(2) unique_roots(2)]*0.14, [0 1],'-')
% plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')
% legend off



%% CD4CD8 type 1 interaction
test=[];
for i=349:363
    test=cat(1,test,CD4(i,312:338,75));
end
test=mean(test,1);

test1=[];
for i=349:363
    test1=cat(1,test1,CD8(i,312:338,75));
end
test1=mean(test1,1);

newcolors = [0 0.7 0.7 ; 0.7 0 0.7 ; 0 0.7 0.7;  1 0 0; 0 0 0 ; 0.7 0 0.7;  1 0 0; 0 0 0  ];
colororder(newcolors)
plot((1:27)*0.14,test/max(test(:)),(1:27)*0.14,test1/max(test1(:)))

[f,gof]=fit((1:27)',test'/max(test(:)),'poly9');
hold on   

xf = 1:27;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:27);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')

[f,gof]=fit((1:27)',test1'/max(test1(:)),'poly9');
hold on   

xf = 1:27;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:27);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')
legend off


%% CD4CD8 type 2 interaction
CD4rotate = imrotate(CD4,35);
CD8rotate = imrotate(CD8,35);
test=[];
for i=456:458
    test=cat(1,test,CD4rotate(i,425:445,78));
end
test=mean(test,1);

test1=[];
for i=456:458
    test1=cat(1,test1,CD8rotate(i,425:445,78));
end
test1=mean(test1,1);

newcolors = [0 0.7 0.7 ; 0.7 0 0.7 ; 0 0.7 0.7;  1 0 0; 0 0 0 ; 0.7 0 0.7;  1 0 0; 0 0 0  ];
colororder(newcolors)
plot((1:21)*0.14,test/max(test(:)),(1:21)*0.14,test1/max(test1(:)))

[f,gof]=fit((1:21)',test'/max(test(:)),'poly9');
hold on   

xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(5)*0.14,f(unique_roots(5)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(5) unique_roots(5)]*0.14, [0 1],'-')

[f,gof]=fit((1:21)',test1'/max(test1(:)),'poly9');
hold on   

xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')
legend off

%% CD4CD8 another type 1 interaction
CD4rotate = imrotate(CD4,85);
CD8rotate = imrotate(CD8,85);
test=[];
for i=367:382
    test=cat(1,test,CD4rotate(i,356:376,82));
end
test=mean(test,1);

test1=[];
for i=367:382
    test1=cat(1,test1,CD8rotate(i,356:376,82));
end
test1=mean(test1,1);

newcolors = [0 0.7 0.7 ; 0.7 0 0.7 ; 0 0.7 0.7;  1 0 0; 0 0 0 ; 0.7 0 0.7;  1 0 0; 0 0 0  ];
colororder(newcolors)
plot((1:21)*0.14,test/max(test(:)),(1:21)*0.14,test1/max(test1(:)))

[f,gof]=fit((1:21)',test'/max(test(:)),'poly9');
hold on   

xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(2)*0.14,f(unique_roots(2)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(2) unique_roots(2)]*0.14, [0 1],'-')

[f,gof]=fit((1:21)',test1'/max(test1(:)),'poly9');
hold on   

xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')
legend off



%% CD8 loose junction
CD8permute = permute(CD8,[3 2 1]);
CD8rotate = imrotate(CD8permute,45);

test=[];
for i=234:238
    test=cat(1,test,CD8rotate(i,326:380,353));
end
test=mean(test,1);


newcolors = [0 0.7 0.7 ; 0 0.8 0 ; 1 0 1];
colororder(newcolors)
plot((1:55)*0.14,test/max(test(:)))


%%%%%%%%%%%%%% - IM Direct interactions activate T cell on tumor cell in IM
MART1 = extractChannels(1,bfGetReader('D:\IM_activatedTcell_tumor_type1.tif'),0,[]);
CD8 = extractChannels(5,bfGetReader('D:\IM_activatedTcell_tumor_type1.tif'),0,[]);

test=[];
for i=461:466
    test=cat(1,test,MART1(i,512:532,97));
end
test=mean(test,1);

test1=[];
for i=461:466
    test1=cat(1,test1,CD8(i,512:532,97));
end
test1=mean(test1,1);
newcolors = [0 1 0 ; 1 0 1 ; 0 1 0 ;  1 0 0; 0 0 0; 1 0 1 ; 1 0 0; 0 0 0];
colororder(newcolors)

[f,gof]=fit((1:21)',test(1:21)'/max(test(:)),'poly9');
colororder(newcolors)
plot((1:21)*0.14, test(1:21)/max(test(:)), (1:21)*0.14, test1(1:21)/max(test1(:))) 
hold on   
xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(3)*0.14,f(unique_roots(3)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(3) unique_roots(3)]*0.14, [0 1],'-')



[f,gof]=fit((1:21)',test1(1:21)'/max(test1(:)),'poly9');
colororder(newcolors)
hold on   
xf = 1:21;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:21);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')

legend off

%% low vs high SNR comparison
PDL1 = extractChannels(16,bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\Dataset1-LSP13626-melanoma_in-situ.tif'),0,[1742 3778 200 200]);
lowPD1 = extractChannels(18,bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\Dataset1-LSP13626-melanoma_in-situ.tif'),0,[1742 3778 200 200]);
highPD1 = extractChannels(40,bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\Dataset1-LSP13626-melanoma_in-situ.tif'),0,[1742 3778 200 200]);
MHC1 = extractChannels(7,bfGetReader('N:\lsp-analysis\cycif-techdev\ClarenceLSM980\20micronF8II\Dataset1-LSP13626-melanoma_in-situ.tif'),0,[1742 3778 200 200]);

PDL1rotate = imrotate(PDL1,45);
highPD1rotate = imrotate(highPD1,45);
lowPD1rotate = imrotate(lowPD1,45);

test=[];
for i=112:135
    test=cat(2,test,PDL1rotate(135:140,i,88));
end
test=mean(test,1);

test1=[];
for i=112:135
    test1=cat(2,test1,highPD1rotate(135:140,i,88));
end
test1=mean(test1,1);

test2=[];
for i=112:135
    test2=cat(2,test2,lowPD1rotate(135:140,i,88));
end
test2=mean(test2,1);
newcolors = [0.8 0.8 0 ; 1 0 1 ; 0 1 0 ; 0.7 0.7 0 ; 1 0 0 ; 0 0 0 ; 1 0 1 ; 1 0 0 ; 0 0 0 ; 0 1 0 ; 1 0 0 ; 0 0 0];
colororder(newcolors)

pixelSize = 0.14;
[f,gof]=fit((1:23)',test(1:23)'/max(test(:)),'poly9');
colororder(newcolors)
plot((1:23)*0.14, test(1:23)/max(test(:)), (1:23)*0.14, test1(1:23)/max(test1(:)), (1:23)*0.14, test2(1:23)/max(test2(:))) 
hold on   
xf = 1:23;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:23);
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
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')
disp(num2str(unique_roots(idx)*pixelSize))

[f,gof]=fit((1:23)',test1(1:23)'/max(test1(:)),'poly9');
colororder(newcolors)
hold on   
xf = 1:23;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:23);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test1)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test1)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(2)*0.14,f(unique_roots(2)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(2) unique_roots(2)]*0.14, [0 1],'-')
maxF=0;
for i=1:numel(unique_roots)
    if f(unique_roots(i))>maxF
        maxF = f(unique_roots(i));
        idx = i;
    end
end
disp(num2str(unique_roots(idx)*pixelSize))

[f,gof]=fit((1:23)',test2(1:23)'/max(test2(:)),'poly9');
colororder(newcolors)
hold on   
xf = 1:23;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')
 [fx fxx]= differentiate(f,1:0.5:23);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(test2)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test2)
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
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')
disp(num2str(unique_roots(idx)*pixelSize))
legend('PDL1','high SNR PD1', 'low SNR PD1')
imshow([imadjust(PDL1(:,:,88),[0 0.05]) imadjust(highPD1(:,:,88),[0 0.2]) imadjust(lowPD1(:,:,88),[0 0.01])])

%% panCK low vs high SNR vs CD11c
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-03_processed.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))

  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-03_processed.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(2000:2380,1420:1860,:,:);
   Irotated = imrotate(I,-20);
   panCK=Irotated(:,:,:,2);
   CD11c=Irotated(:,:,:,5);
   lowpanCK = Irotated(:,:,:,4);
    
scalingFactor=1;
pixelSize=0.06/scalingFactor;
xlow = round(142*scalingFactor);
xhigh = round(180*scalingFactor);
ylow = round(199*scalingFactor);
yhigh = round(212*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,panCK(ylow:yhigh,i,45));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD11c(ylow:yhigh,i,45));
end
test1=mean(test1,1);

test2=[];
for i=xlow:xhigh
    test2=cat(2,test2,lowpanCK(ylow:yhigh,i,45));
end
test2=mean(test2,1);


newcolors = [0.7 0.7 0 ; 1 0 1 ; 0 1 0 ; 0.7 0.7 0 ; 1 0 0 ; 0 0 0 ; 1 0 1 ; 1 0 0 ; 0 0 0 ; 0 1 0 ; 1 0 0 ; 0 0 0];
colororder(newcolors)


[f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
colororder(newcolors)
plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:)), (1:length)*pixelSize, test2(1:length)/max(test2(:))) 
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
disp(num2str(unique_roots(idx)*pixelSize))

[f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
colororder(newcolors)
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
disp(num2str(unique_roots(idx)*pixelSize))

[f,gof]=fit((1:length)',test2(1:length)'/max(test2(:)),'poly9');
colororder(newcolors)
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
              for i = 1:numel(test2)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test2)
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
disp(num2str(unique_roots(idx)*pixelSize))
%% low vs high resolution on non-interacting CD3 CD11c
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-03_processed.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))

  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-03_processed.ome.tiff',1);
   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end
   Icrop=vol(1380:2000,1640:1970,:,:);
   Irotated = imrotate(Icrop,-45);


scalingFactor=0.05;
pixelSize=0.06/scalingFactor;
CD3=Irotated(:,:,:,3);
CD11c=Irotated(:,:,:,5);
CD3 = imresize(CD3,scalingFactor,'nearest');
CD11c = imresize(CD11c,scalingFactor,'nearest');

xlow = round(425*scalingFactor);
xhigh = round(525*scalingFactor);
ylow = round(216*scalingFactor);
yhigh = round(232*scalingFactor);
length = xhigh-xlow+1;

test=[];
for i=xlow:xhigh
    test=cat(2,test,CD3(ylow:yhigh,i,55));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD11c(ylow:yhigh,i,55));
end
test1=mean(test1,1);

newcolors = [0 1 0 ; 1 0 1 ; 0 1 0 ;  1 0 0; 0 0 0; 1 0 1 ; 1 0 0; 0 0 0];
colororder(newcolors)



[f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly7');
colororder(newcolors)
plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
disp(num2str(unique_roots(idx)*pixelSize))

[f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly7');
colororder(newcolors)
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
disp(num2str(unique_roots(idx)*pixelSize))


%% low vs high resolution on interacting panCK vs CD11c
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-03_processed.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))

  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-03_processed.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(2000:2380,1420:1860,:,:);
   Irotated = imrotate(I,-20);


    
scalingFactor=1;
panCK=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
CD11c=imresize(Irotated(:,:,:,5),scalingFactor,'nearest');
pixelSize=0.06/scalingFactor;
xlow = round(130*scalingFactor);
xhigh = round(195*scalingFactor);
ylow = round(199*scalingFactor);
yhigh = round(212*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,panCK(ylow:yhigh,i,45));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD11c(ylow:yhigh,i,45));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly5');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly5');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))
legend('panCK','CD11c')

%% panCK low vs high SNR vs CD11c
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-04_processed-Stitching-01_2.tif');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))

Hoechst = extractChannelsChannelZ(1,bfGetReader('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-04_processed-Stitching-01_2.ome.tif'),0,[]);
panCK = extractChannelsChannelZ(2,bfGetReader('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-04_processed-Stitching-01_2.tif'),0,[]);
lowpanCK = extractChannelsChannelZ(4,bfGetReader('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-04_processed-Stitching-01_2.tif'),0,[]);
CD11c = extractChannelsChannelZ(5,bfGetReader('U:\cycif-techdev\Clarence_LSM980\interactions\panCKCD103CD11c\cycle1-04_processed-Stitching-01_2.tif'),0,[]);
    
scalingFactor=1;
pixelSize=0.06/scalingFactor;
xlow = round(347*scalingFactor);
xhigh = round(380*scalingFactor);
ylow = round(276*scalingFactor);
yhigh = round(276*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,panCK(ylow:yhigh,i,30));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD11c(ylow:yhigh,i,30));
end
test1=mean(test1,1);

test2=[];
for i=xlow:xhigh
    test2=cat(2,test2,lowpanCK(ylow:yhigh,i,30));
end
test2=mean(test2,1);


newcolors = [0.7 0.7 0 ; 1 0 1 ; 0 1 0 ; 0.7 0.7 0 ; 1 0 0 ; 0 0 0 ; 1 0 1 ; 1 0 0 ; 0 0 0 ; 0 1 0 ; 1 0 0 ; 0 0 0];
colororder(newcolors)


[f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
colororder(newcolors)
plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:)), (1:length)*pixelSize, test2(1:length)/max(test2(:))) 
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
disp(num2str(unique_roots(idx)*pixelSize))

[f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
colororder(newcolors)
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
disp(num2str(unique_roots(idx)*pixelSize))

[f,gof]=fit((1:length)',test2(1:length)'/max(test2(:)),'poly9');
colororder(newcolors)
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
              for i = 1:numel(test2)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(test2)
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
disp(num2str(unique_roots(idx)*pixelSize))
   

%% low vs high resolution on interacting panCK vs CD11c
%40x confocal
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-18\New-01_processed.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))

  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-18\New-01_processed.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(953:1271,908:1228,:,:);
   Irotated = imrotate(I,-35);
    
scalingFactor=1;
CD4=imresize(Irotated(:,:,:,1),scalingFactor,'nearest');
CD8=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
pixelSize=0.08/scalingFactor;
xlow = round(175*scalingFactor);
xhigh = round(225*scalingFactor);
ylow = round(169*scalingFactor);
yhigh = round(169*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD4(ylow:yhigh,i,35));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD8(ylow:yhigh,i,35));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))


%20x widefield
     metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-18\20x WF-06.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))

  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-18\20x WF-06.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(992:1113,908:1125,:,:);
   Irotated = imrotate(I,-25);
    
scalingFactor=0.65;
CD4=imresize(Irotated(:,:,:,1),scalingFactor,'nearest');
CD8=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
pixelSize=0.21/scalingFactor;
xlow = round(138*scalingFactor);
xhigh = round(167*scalingFactor);
ylow = round(94*scalingFactor);
yhigh = round(94*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD4(ylow:yhigh,i,6));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD8(ylow:yhigh,i,6));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly7');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly7');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))

%% low vs high resolution on interacting CD4 vs CD8
%40x confocal
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\40x confocal-01_processed.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))

  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\40x confocal-01_processed.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(1080:1422,1449:1836,:,:);
   Irotated = imrotate(I,-35);
    
scalingFactor=1;
CD4=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
CD8=imresize(Irotated(:,:,:,3),scalingFactor,'nearest');
pixelSize=0.08/scalingFactor;
xlow = round(258*scalingFactor);
xhigh = round(305*scalingFactor);
ylow = round(228*scalingFactor);
yhigh = round(232*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD4(ylow:yhigh,i,24));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD8(ylow:yhigh,i,24));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))
    legend('CD4','CD8')
imshowpair(imadjust(I(:,:,24,2),[0.01 0.3]),imadjust(I(:,:,24,3),[0.01 0.07]))
%20x widefield
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\20x widefield.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))
  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\20x widefield.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(1393:1567,930:1158,:,:);
   Irotated = imrotate(I,-35);
    
scalingFactor=1;
CD4=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
CD8=imresize(Irotated(:,:,:,3),scalingFactor,'nearest');
pixelSize=0.170/scalingFactor;
xlow = round(150*scalingFactor);
xhigh = round(185*scalingFactor);
ylow = round(123*scalingFactor);
yhigh = round(127*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD4(ylow:yhigh,i,9));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD8(ylow:yhigh,i,9));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))


    %20x confocal
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\20x confocal-02-LSM Plus Processing-03.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))
  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\20x confocal-02-LSM Plus Processing-03.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(1261:1397,815:971,:,:);
   Irotated = imrotate(I,-35);
    
scalingFactor=1;
CD4=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
CD8=imresize(Irotated(:,:,:,3),scalingFactor,'nearest');
pixelSize=0.210/scalingFactor;
xlow = round(92*scalingFactor);
xhigh = round(116*scalingFactor);
ylow = round(92*scalingFactor);
yhigh = round(97*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD4(ylow:yhigh,i,9));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD8(ylow:yhigh,i,9));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))
legend('CD4','CD8')
imshowpair(imadjust(I(:,:,9,2),[0.01 0.3]),imadjust(I(:,:,9,3),[0.01 0.07]))

      %40x widefield
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\separate\40x widefield-01.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))
  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\separate\40x widefield-01.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(781:1113,1063:1503,:,:);
   Irotated = imrotate(I,-35);
    
scalingFactor=1;
CD4=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
CD8=imresize(Irotated(:,:,:,3),scalingFactor,'nearest');
pixelSize=0.086/scalingFactor;
xlow = round(289*scalingFactor);
xhigh = round(335*scalingFactor);
ylow = round(230*scalingFactor);
yhigh = round(230*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD4(ylow:yhigh,i,10));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD8(ylow:yhigh,i,10));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))
    legend('CD4','CD8')
    imshowpair(imadjust(I(:,:,10,2),[0 0.05]),imadjust(I(:,:,10,3),[0 0.02]))

       %% single vs multiple lines 40x confocal
 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\separate\40x confocal-03_processed.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))
  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\separate\40x confocal-03_processed.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(1078:1381,1341:1710,:,:);
   Irotated = imrotate(I,35);
    
scalingFactor=1;
CD4=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
CD8=imresize(Irotated(:,:,:,3),scalingFactor,'nearest');
pixelSize=0.08/scalingFactor;
xlow = round(173*scalingFactor);
xhigh = round(231*scalingFactor);
ylow = round(270*scalingFactor);
yhigh = round(270*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD4(ylow:yhigh,i,10));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,CD8(ylow:yhigh,i,10));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))
legend('CD4','CD8')
imshowpair(imadjust(I(:,:,10,2),[0 0.2]),imadjust(I(:,:,10,3),[0 0.08]))


%% multiple z planes in 40x confocal

 metadata = loci.formats.Memoizer(bfGetReader());
    options = javaObject('loci.formats.in.DynamicMetadataOptions');
  
    metadata.setMetadataOptions(options);
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\40x confocal-01_processed.ome.tiff');
    sizeX = metadata.getSizeX;
    sizeY = metadata.getSizeY;
    sizeZ = metadata.getSizeZ;
    numChan = metadata.getSizeC;
disp(num2str(metadata.getSeriesCount()))

  [I] = bftileopen('U:\cycif-techdev\Clarence_LSM980\CD4CD8interactions\2024-08-19\40x confocal-01_processed.ome.tiff',1);

   vol = zeros(sizeY,sizeX,sizeZ,numChan,'uint16');
   for iC = 1:numChan
       for iZ =  1:sizeZ
            vol(:,:,iZ,iC) = I{1,1}{(iC-1)*sizeZ+iZ};
       end
   end

   I=vol(1080:1422,1449:1836,:,:);
   Irotated = imrotate(I,-35);
    
scalingFactor=1;
CD4=imresize(Irotated(:,:,:,2),scalingFactor,'nearest');
CD8=imresize(Irotated(:,:,:,3),scalingFactor,'nearest');
pixelSize=0.08/scalingFactor;
xlow = round(258*scalingFactor);
xhigh = round(305*scalingFactor);
ylow = round(230*scalingFactor);
yhigh = round(230*scalingFactor);
length = xhigh-xlow+1;
test=[];
for iZ = 24:24
    test=cat(3,test,CD4(ylow:yhigh,xlow:xhigh,iZ));
end
test=mean(test,3);

test1=[];
for iZ = 24:24
    test1=cat(3,test1,CD8(ylow:yhigh,xlow:xhigh,iZ));
end
test1=mean(test1,3);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))
legend('CD4','CD8')
imshowpair(imadjust(max(I(:,:,24:24,2),[],3),[0 0.2]),imadjust(max(I(:,:,24:24,3),[],3),[0 0.04]))


%% CHTN cell interaction 1
metadata = loci.formats.Memoizer(bfGetReader());
    metadata.setId('U:\cycif-techdev\Clarence_LSM980\20umCHTN\CHTN16bitbgsub.tif')  

coords = [245 1300 295 450; 880 1421 1042 691];
cells = [];    
for iCell = 1:4
    channel = [];
    for iChan = [2]%[3 4 14 15 17 21]
        channel= cat(4,channel,extractChannels(iChan,metadata,0,round([coords(1,iCell) coords(2,iCell) 50 50]/0.14)));
    end
cells=cat(5,cells,channel);
end


metadata = loci.formats.Memoizer(bfGetReader());
metadata.setId('U:\cycif-techdev\Clarence_LSM980\20umCHTN\CHTN16bitbgsub.tif')  
CD3=extractChannels(3,metadata,0,round([245 880 28 28]/0.14));
MART1=extractChannels(4,metadata,0,round([245 880 28 28]/0.14));
CD3rotated=imrotate(CD3,-90);
MART1rotated=imrotate(MART1,-90);

scalingFactor=1;
pixelSize=0.14/scalingFactor;
xlow = round(69*scalingFactor);
xhigh = round(95*scalingFactor);
ylow = round(148*scalingFactor);
yhigh = round(159*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD3rotated(ylow:yhigh,i,54));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,MART1rotated(ylow:yhigh,i,54));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))
    legend('CD3','MART1')
    


%% CHTN cell interaction 4


metadata = loci.formats.Memoizer(bfGetReader());
metadata.setId('U:\cycif-techdev\Clarence_LSM980\20umCHTN\CHTN16bitbgsub.tif')  
CD3=extractChannels(3,metadata,0,round([450 690 28 28]/0.14));
MART1=extractChannels(4,metadata,0,round([450 690 28 28]/0.14));
CD3rotated=imrotate(CD3,-25);
MART1rotated=imrotate(MART1,-25);

scalingFactor=1;
pixelSize=0.14/scalingFactor;
xlow = round(121*scalingFactor);
xhigh = round(134*scalingFactor);
ylow = round(109*scalingFactor);
yhigh = round(112*scalingFactor);
length = xhigh-xlow+1;
test=[];
for i=xlow:xhigh
    test=cat(2,test,CD3rotated(ylow:yhigh,i,54));
end
test=mean(test,1);

test1=[];
for i=xlow:xhigh
    test1=cat(2,test1,MART1rotated(ylow:yhigh,i,54));
end
test1=mean(test1,1);


newcolors = [0 1 0 ; 1 0 1; 0 1 0; 0 0 0; 0 0 0; 1 0 1 ; 0 0 0; 0 0 0];
    colororder(newcolors)
    
    
    [f,gof]=fit((1:length)',test(1:length)'/max(test(:)),'poly9');
    colororder(newcolors)
    plot((1:length)*pixelSize, test(1:length)/max(test(:)), (1:length)*pixelSize, test1(1:length)/max(test1(:))) 
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
            maxF = f(unique_roots(i))
            idx = i;
        end
    end
    plot(unique_roots(idx)*pixelSize,f(unique_roots(idx)),'x','MarkerEdgeColor','red','MarkerSize',12)
    plot([unique_roots(idx) unique_roots(idx)]*pixelSize, [0 1],'-')
        disp(num2str(unique_roots(idx)*pixelSize))
    
    [f,gof]=fit((1:length)',test1(1:length)'/max(test1(:)),'poly9');
    colororder(newcolors)
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
    disp(num2str(unique_roots(idx)*pixelSize))
    legend('CD3','MART1')
    











