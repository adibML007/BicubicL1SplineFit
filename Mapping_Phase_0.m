function [i_bar,j_bar]=Mapping_Phase_0(x,y,xM,yM,I,J)
if xM==max(x)&&yM==max(y)
    i_bar=I;
    j_bar=J;
elseif xM==max(x)&&yM~=max(y)
    i_bar=I;
    [~,j_bar]=Mapping_Phase_1(x,y,xM,yM);
elseif xM~=max(x)&&yM==max(y)
    j_bar=I;
    [i_bar,~]=Mapping_Phase_1(x,y,xM,yM);
else
    [i_bar,j_bar]=Mapping_Phase_1(x,y,xM,yM);
end