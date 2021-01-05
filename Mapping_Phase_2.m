function [i,j,n]=Mapping_Phase_2(i_bar,j_bar,I,J)

if i_bar==I&&j_bar==J           %Last node
    n=(I-1)*(J-1);
    i=I;
    j=J;
elseif i_bar==I&&j_bar~=J       %Boundary nodes in x direction
    n=(J-1)*(i_bar-2)+j_bar;
    i=I;
    j=j_bar;
elseif i_bar~=I&&j_bar==J       %Boundary nodes in y direction
    n=i_bar*(J-1);
    i=i_bar;
    j=J;
else
    n=(i_bar-1)*(J-1)+j_bar;    %All other nodes
    i=i_bar;
    j=j_bar;
    
end
end