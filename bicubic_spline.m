function[zplot,cpu_time]=bicubic_spline(z_given)
%This function implements the bicubic spline method to interpolate data
%points based on a given set
%% Input Parameters
%z_given denotes the given data points
%% Output Parameters
%zplot is the interpolated points
tic
x=0:size(z_given,1)-1;
y=0:size(z_given,2)-1;
z=z_given;
c=Bicubic_interpolation(x,y,z);
I=find(x==max(x));
J=find(y==max(y));
zplot=zeros();
%Finding interpolating z-values
xplot=linspace(0,x(end),200);
yplot=linspace(0,y(end),200);
for i=1:numel(xplot)
    for j=1:numel(yplot)
        [i_temp0,j_temp0]=Mapping_Phase_1(x,y,xplot(i),yplot(j));
        [~,p1]=bicubic_grid((xplot(i)-x(i_temp0))/(x(i_temp0+1)-x(i_temp0))...
            ,(yplot(j)-y(j_temp0))/(y(j_temp0+1)-y(j_temp0)));
        [~,~,n]=Mapping_Phase_2(i_temp0,j_temp0,I,J);
        p_n=p1(1,:)';
        c_n=c(1:16,n)';
        zplot_n=c_n*p_n;
        zplot(i,j)=zplot_n;
    end
end
toc
cpu_time=toc;
end