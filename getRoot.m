function idx=getRoot(x,f)

maxF=[];
for i=1:numel(x)
    if f(x(i))>maxF
        maxF = f(x(i));
        idx = i;
    end
end


end