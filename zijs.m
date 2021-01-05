function [zf,zfx,zfy,zfxfy,Zx,Zy,Zxy]=zijs(x,y,z)
zf=z;
zfx=[];
zfy=[];
zfxfy=[];
Zx=[];
Zy=[];
Zxy=[];
for i=1:length(y)  %Number of columns in zf
    [G,X,~]=Spline_interpolation_2nd_derivative(x',zf(:,i));
    zfx=cat(2,zfx,G');
    Zx=cat(1,Zx,X);
end


for j=1:length(x)   %Number of rows in zf
    [F,Y,~]=Spline_interpolation_2nd_derivative(y,zf(j,:));
    zfy=cat(1,zfy,F);
    Zy=cat(1,Zy,Y);
end


for i=1:length(y)
    [GF,XY,~]=Spline_interpolation_2nd_derivative(x',zfy(:,i));
    zfxfy=cat(2,zfxfy,GF');
    Zxy=cat(1,Zxy,XY*diag(Zy((0:length(y):length(y)*(length(x)-1))+i,i)));
end
% disp(Zx)
% disp(XY)
% disp(Zy)
% disp(Zxy)