function [xy,x1]=bicubic_grid(x,y)
%%
%Initialization of different matrices
xy_f=zeros();
xy_fx=zeros();
xy_fy=zeros();
xy_fxy=zeros();

%% Calculation of bicubic grid functions
%Calculating the function values The problem here is that 0^(-ve) is taken
%as Inf in MATLAB which is theoretically 0 in our case. So,we need to
%handle those cases differently when either x or y is zero by manually
%inputting zeros while this type of case occurs.

for s=1:4
    for t=1:4
        if s==1&&t~=1
            xy_f(s,t)=x^(s-1)*y^(t-1);
            xy_fx(s,t)=0;                       %0^(-ve) appears
            xy_fy(s,t)=(t-1)*x^(s-1)*y^(t-2);
            xy_fxy(s,t)=0;                      %0^(-ve) appears
        elseif s~=1&&t==1
            xy_f(s,t)=x^(s-1)*y^(t-1);
            xy_fx(s,t)=(s-1)*x^(s-2)*y^(t-1);
            xy_fy(s,t)=0;                       %0^(-ve) appears
            xy_fxy(s,t)=0;                      %0^(-ve) appears
        elseif s==1&&t==1
            xy_f(s,t)=x^(s-1)*y^(t-1);
            xy_fx(s,t)=0;                       %0^(-ve) appears
            xy_fy(s,t)=0;                       %0^(-ve) appears
            xy_fxy(s,t)=0;                      %0^(-ve) appears
        else
            xy_f(s,t)=x^(s-1)*y^(t-1);
            xy_fx(s,t)=(s-1)*x^(s-2)*y^(t-1);
            xy_fy(s,t)=(t-1)*x^(s-1)*y^(t-2);
            xy_fxy(s,t)=(s-1)*(t-1)*x^(s-2)*y^(t-2);
        end
    end
end


%Construction of the final matrix involving f,fx,fy and fxy for each node
xy=cat(1,reshape(xy_f',1,length(xy_f)^2),reshape(xy_fx',1,length(xy_fx)^2),...
    reshape(xy_fy',1,length(xy_fy)^2),reshape(xy_fxy',1,length(xy_fxy)^2));
x1=xy(1,:);
end
