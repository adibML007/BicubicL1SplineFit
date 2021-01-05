function [X,Y,Z]=drawspline2D(x,y,z,px,qy)
c=Bicubic_interpolation(x,y,z);
I=find(x==max(x));
J=find(y==max(y));
zplot=zeros();
%prompt=input('Input 1 to display given and interpolated node values to check for consistency\n');
% prompt=1;
% if prompt==1
%     q1=2;
%     q2=2;
% else
if nargin==3
    q1=input('Please, input the number of interpolation points for x direction\n');
    q2=input('Please, input the number of interpolation points for y direction\n');
else
    q1=px;
    q2=qy;
end
% end

%% Parameters for smoothing splines
a1=min(x);
n1=numel(x)-1;
h1=(max(x)-min(x))/n1;
a2=min(y);
n2=numel(y)-1;
h2=(max(y)-min(y))/n2;
s10=zeros();
s20=zeros();

%% Smooth 2D spline surface plot
%Creating matrices of data points in both direction
for w1=1:n1
    x1=linspace((a1+(w1-1)*h1),(a1+w1*h1),q1);
    s10(w1,1:numel(x1))=x1;
end

for w2=1:n2
    x2=linspace((a2+(w2-1)*h2),(a2+w2*h2),q2);
    s20(w2,1:numel(x2))=x2;
end

%Reshaping matrices into vectors
s11=reshape(s10',numel(s10),1);
s22=reshape(s20',numel(s20),1);
xplot=unique(s11);
yplot=unique(s22);

%Finding interpolating z-values
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
% if q1==2&&q2==2
%     disp(z)
%     disp(zplot)
%     sum(sum(z==zplot))
% else
if nargout>3
    disp('Now, the interpolated 2-D spline surfaces will be plotted.')
    %Plotting Parameters
    if min(x)<=min(xplot)
        Xmin=min(x);
    else
        Xmin=min(xplot);
    end
    if max(x)>=max(xplot)
        Xmax=max(x);
    else
        Xmax=max(xplot);
    end
    if min(y)<=min(yplot)
        Ymin=min(y);
    else
        Ymin=min(yplot);
    end
    if max(y)>=max(yplot)
        Ymax=max(y);
    else
        Ymax=max(yplot);
    end
    if min(min(z))<=min(min(zplot))
        Zmin=min(min(z));
    else
        Zmin=min(min(zplot));
    end
    if max(max(z))>=max(max(zplot))
        Zmax=max(max(z));
    else
        Zmax=max(max(zplot));
    end
    ax1=subplot(1,2,1);
    surf(x,y,z')
    title('Given function value curve')
    ax2=subplot(1,2,2);
    surf(xplot,yplot,zplot')
    title('Approximated function value curve')
    axis([ax1 ax2],[Xmin Xmax Ymin Ymax Zmin Zmax])
    shading interp
else
    X=xplot;
    Y=yplot;
    Z=zplot;
end