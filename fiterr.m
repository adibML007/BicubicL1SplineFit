function[e]=fiterr(xM,yM,zM,xMyM,x,y,I,J,c)
error=zeros();
for m=1:length(xM)
    [i_bar,j_bar]=Mapping_Phase_0(x,y,xM(m),yM(m),I,J);
    [~,~,n]=Mapping_Phase_2(i_bar,j_bar,I,J);
    p_n=xMyM(m,1:16)';
    c_n=c(1:16,n)';
    error(m)=abs(c_n*p_n-zM(m));
end
e=sum(error);
end