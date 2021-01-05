function [i,j]=Mapping_Phase_1(x,y,xM,yM)

for k=1:(length(x)-1)
    if xM>=x(k)&&xM<=x(k+1)
        i=k;
    end
end

for k=1:(length(y)-1)
    if yM>=y(k)&&yM<=y(k+1)
        j=k;
    end
end
end