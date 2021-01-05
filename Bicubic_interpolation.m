function [c,Zx,Zy,Zxy]=Bicubic_interpolation(x,y,z0)
c=zeros();
[zf,zfx,zfy,zfxfy,Zx,Zy,Zxy]=zijs(x,y,z0);
n=1;
%Define the constant coefficient 16 by 16 matrix p
p=[  1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0
     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0
     1     0     0     0     1     0     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     1     0     0     0     2     0     0     0     3     0     0     0
     0     1     0     0     0     1     0     0     0     1     0     0     0     1     0     0
     0     0     0     0     0     1     0     0     0     2     0     0     0     3     0     0
     1     1     1     1     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     1     1     1     1     0     0     0     0     0     0     0     0
     0     1     2     3     0     0     0     0     0     0     0     0     0     0     0     0
     0     0     0     0     0     1     2     3     0     0     0     0     0     0     0     0
     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
     0     0     0     0     1     1     1     1     2     2     2     2     3     3     3     3
     0     1     2     3     0     1     2     3     0     1     2     3     0     1     2     3
     0     0     0     0     0     1     2     3     0     2     4     6     0     3     6     9];

for i=1:(length(x)-1)
    for j=1:(length(y)-1)
        z_p=[zf(i,j);zfx(i,j);zfy(i,j);zfxfy(i,j);zf(i+1,j);zfx(i+1,j);...
            zfy(i+1,j);zfxfy(i+1,j);zf(i,j+1);zfx(i,j+1);zfy(i,j+1);...
            zfxfy(i,j+1);zf(i+1,j+1);zfx(i+1,j+1);zfy(i+1,j+1);...
            zfxfy(i+1,j+1)];
        c(1:16,n)=p\z_p;
        n=n+1;
    end
end
end