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


%% MART1 CD4 type 2 interaction
MART1permute = permute(MART1,[3 2 1]);
CD4permute = permute(CD4,[3 2 1]);
test=[];
for i=68:72
    test=cat(1,test,CD4permute(i,257:287,352));
end
test=mean(test,1);

test1=[];
for i=68:72
    test1=cat(1,test1,MART1permute(i,257:287,352));
end
test1=mean(test1,1);

newcolors = [0 0.7 0.7 ; 0 0.8 0 ; 1 0 1];
colororder(newcolors)
plot((1:31)*0.14,test/max(test(:)),(1:31)*0.14,test1/max(test1(:)))


lineintegral = max(cat(1,test/max(test(:)), test1/max(test1(:))));
[f,gof]=fit((1:31)',lineintegral','poly9');
newcolors = [0 0.7 0.7 ; 0 0.8 0 ; 0 0 0; 0 0 0; 0 0 0;0 0 0; 0 0 0];
colororder(newcolors)
plot((1:28)*0.14, test(1:28)/max(test(:)), (1:28)*0.14, test1(1:28)/max(test1(:))) 
hold on   

xf = 1:28;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:31);
              fdxHan=@(x)(differentiate(f,x));
 dxroots =[];
              for i = 1:numel(lineintegral)
                dxroots(i)=fzero(fdxHan, i);
                if dxroots(i) <0 || dxroots(i) > numel(lineintegral)
                    dxroots(i)=NaN;
                end
              end
unique_roots=uniquetol(dxroots,0.001);
plot(unique_roots(2)*0.14,f(unique_roots(2)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot(unique_roots(4)*0.14,f(unique_roots(4)),'x','MarkerEdgeColor','red','MarkerSize',12)
plot([unique_roots(2) unique_roots(2)]*0.14, [0 1],'-')
plot([unique_roots(4) unique_roots(4)]*0.14, [0 1],'-')
legend off



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
for i=458:463
    test=cat(1,test,CD4rotate(i,420:450,78));
end
test=mean(test,1);

test1=[];
for i=458:463
    test1=cat(1,test1,CD8rotate(i,420:450,78));
end
test1=mean(test1,1);

newcolors = [0 0.7 0.7 ; 0.7 0 0.7 ; 0 0.7 0.7;  1 0 0; 0 0 0 ; 0.7 0 0.7;  1 0 0; 0 0 0  ];
colororder(newcolors)
plot((1:31)*0.14,test/max(test(:)),(1:31)*0.14,test1/max(test1(:)))

[f,gof]=fit((1:31)',test'/max(test(:)),'poly9');
hold on   

xf = 1:31;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:31);
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

[f,gof]=fit((1:31)',test1'/max(test1(:)),'poly9');
hold on   

xf = 1:31;
yf=[];
for i=xf
    yf=cat(2,yf,f(xf(i)));
end
plot(xf*0.14,yf,'--')

 [fx fxx]= differentiate(f,1:0.5:31);
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
legend off




%%CD8 loose junction
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


