clear
xM=[-2 -1.5 1.2 3.5 5];
yM=xM;
x=xM;
y=yM;
z=zeros();
X=1:11;
disp('Please see the following table for different function options\n')
F={'Linear','Polynomial','Sinusoidal','Rosenbrock','Absolute','Random',...
    'Quotient function','Baltimore','Killeen','Peaks','End Program'};
t=table(X',F','VariableNames',{'Option' 'Function'});
disp(t)
sel=input('Enter your choice: ');
switch sel
    case 1
        z0=@(m,n) m+n;                      %Linear
    case 2
        z0=@(m,n) 0.14*m^4+0.2*m^3*n^2-0.15*m*n^2+0.3*n^3-13;%Polynomial 
    case 3
        z0=@(m,n) (1-4*m-m.^3/17)*sin(n^2); %Sinusoidal 
    case 4
        z0=@(m,n) (1-m)^2+100*(n-m^2)^2;    %Rosenbrock 
    case 5
        z0=@(m,n) abs(m-n);                 %Absolute 
    case 6
        z0=@(m,n) rand();                   %Random 
    case 7
        z0=@(m,n) (m+n)/(m^2+n^2+1);        %Quotient
    case 8
        %Given subset of Baltimore_1000x1000 data
        load('C:\Users\Adib\Downloads\Baltimore_1000x1000.mat')
        x1=510;
        x2=595;
        y1=460;
        y2=500;
        xM=(0:length(x1:2:x2)-1);
        yM=(0:length(y1:2:y2)-1);
        x=xM;
        y=yM;
        z0=G(x1:2:x2,y1:2:y2,:);
        z=z0;
    case 9
        %Given subset of Killeen_901x901 data
%         load('C:\Users\Adib\Downloads\Killeen_901x901.mat')
%         load('E:\01-14-2018\self.mat')
        x1=1;
        x2=600;
        y1=1;
        y2=600;
        xM=(0:length(x1:x2)-1);
        yM=(0:length(y1:y2)-1);
        x=xM;
        y=yM;
        z0=self(x1:x2,y1:y2,:);
        z=z0;
    case 10
        z0=@(m,n) 3*(1-m).^2.*exp(-(m.^2) - (n+1).^2)-...
            10*(m/5 - m.^3 - n.^5).*exp(-m.^2-n.^2)...
            -1/3*exp(-(m+1).^2 - n.^2);
        x=-4:4;
        y=x;
    otherwise 
        z0=@(m,n) m/(n+1);
end

if sel~=8&&sel~=9
for i=1:length(x)
    for j=1:length(y)
            z(i,j)=z0(x(i),y(j));
    end
end
end
z_org=zeros();
[zf,zfx,zfy,zfxfy,~,~,~]=zijs(x,y,z);
c=Bicubic_interpolation(x,y,z);
I=find(x==max(x));
J=find(y==max(y));
xplot=linspace(min(x),max(x),20);
yplot=linspace(min(y),max(y),20);
zplot=zeros();
    for i=1:numel(xplot)
        for j=1:numel(yplot)
            [i_temp0,j_temp0]=Mapping_Phase_1(x,y,xplot(i),yplot(j));
            [~,p1]=bicubic_grid((xplot(i)-x(i_temp0))/(x(i_temp0+1)...
                -x(i_temp0)),(yplot(j)-y(j_temp0))/...
                (y(j_temp0+1)-y(j_temp0)));
            [~,~,n]=Mapping_Phase_2(i_temp0,j_temp0,I,J);
            p_n=p1(1,:)';
            c_n=c(1:16,n)';
            zplot_n=c_n*p_n;
            zplot(i,j)=zplot_n;
        end
    end
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
title('Given Surface') 
ax2=subplot(1,2,2); 
surf(xplot,yplot,zplot')
title('Interpolated Bicubic Spline Surface')
axis([ax1 ax2],[Xmin Xmax Ymin Ymax Zmin Zmax])
% savefile='testResults.mat';
% save(savefile,'c','x','y','z','xplot','yplot','zplot');
% storedStructure = load('testResults.mat');
% imageArray = storedStructure.z;
% figure(1);
% subplot(1,2,1),imshow(imageArray, [])
% imageArray = storedStructure.zplot;
% subplot(1,2,2),imshow(imageArray, [])
