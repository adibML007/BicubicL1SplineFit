%% Function to calculate numerical gradient

%This function calculates numerical gradient for 2-D cubic L1 spline fits
%and store some values in matrix forms to be used later in the steepest
%descent algorithm a.k.a. gradient search

%% Description of function arguments

%Input arguments
%x is the vector for horizontal grid indices
%y is the vector for horizontal grid indices
%z0 are spline node values at each iteration in matrix form
%xM is a vector of given x-coordinate values
%yM is a vector of given y-coordinate values
%zM is a vector of given z values that would be approximated

%Output arguments
%delzij are the gradient values for each node in matrix form
%c is the coefficient matrix for each bicubic interpolation in the grid
%epsilon is the absolute error for each given data point at each iteration

%% Start of the main function
function[delzij]=Gradient_finder(x,y,z0,xM,yM,zM,xMyM)
dgstdzij_temp=[1 0 -3 2 0 0 0 0 -3 0 9 -6 2 0 -6 4]';
dgstdzi1j_temp=[0 0 0 0 0 0 0 0 3 0 -9 6 -2 0 6 -4]';
dgstdzij1_temp=[0 0 3 -2 0 0 0 0 0 0 -9 6 0 0 6 -4]';
dgstdzi1j1_temp=[0 0 0 0 0 0 0 0 0 0 9 -6 0 0 -6 4]';
dgstdzxij_temp=[0 0 0 0 1 0 -3 2 -2 0 6 -4 1 0 -3 2]';
dgstdzxi1j_temp=[0 0 0 0 0 0 0 0 -1 0 3 -2 1 0 -3 2]';
dgstdzxij1_temp=[0 0 0 0 0 0 3 -2 0 0 -6 4 0 0 3 -2]';
dgstdzxi1j1_temp=[0 0 0 0 0 0 0 0 0 0 -3 2 0 0 3 -2]';
dgstdzyij_temp=[0 1 -2 1 0 0 0 0 0 -3 6 -3 0 2 -4 2]';
dgstdzyi1j_temp=[0 0 0 0 0 0 0 0 0 3 -6 3 0 -2 4 -2]';
dgstdzyij1_temp=[0 0 -1 1 0 0 0 0 0 0 3 -3 0 0 -2 2]';
dgstdzyi1j1_temp=[0 0 0 0 0 0 0 0 0 0 -3 3 0 0 2 -2]';
dgstdzxyij_temp=[0 0 0 0 0 1 -2 1 0 -2 4 -2 0 1 -2 1]';
dgstdzxyi1j_temp=[0 0 0 0 0 0 0 0 0 -1 2 -1 0 1 -2 1]';
dgstdzxyij1_temp=[0 0 0 0 0 0 -1 1 0 0 2 -2 0 0 -1 1]';
dgstdzxyi1j1_temp=[0 0 0 0 0 0 0 0 0 0 1 -1 0 0 -1 1]';
[c,Zx,Zy,Zxy]=Bicubic_interpolation(x,y,z0);
I=find(x==max(x));
J=find(y==max(y));
delzij=zeros();
delzij_temp=zeros();
%% Calculate the gradient of each node at a time with all given data points
for i=1:length(x)
    for j=1:length(y)
        for m=1:length(xM)
            %             %% Mapping Phase 1: (xM,yM) to (i,j)
            [i_temp0,j_temp0]=Mapping_Phase_1(x,y,xM(m),yM(m));
            [i_bar,j_bar]=Mapping_Phase_0(x,y,xM(m),yM(m),I,J);
            %% Mapping Phase 2: (i,j) to n
            [~,~,n]=Mapping_Phase_2(i_bar,j_bar,I,J);
            %% Finding the partial derivatives
            dzxdzij_temp1=Zx(1+(j-1)*(length(x)):1+(j-1)*...
                (length(x))+(length(x)-1),1:length(x));
            dzxdzij=[dzxdzij_temp1(i_temp0,i)*(j==j_temp0)...
                dzxdzij_temp1(i_temp0+1,i)*(j==j_temp0)...
                dzxdzij_temp1(i_temp0,i)*(j==j_temp0+1)...
                dzxdzij_temp1(i_temp0+1,i)*(j==j_temp0+1)];
            
            dzydzij_temp1=Zy(1+(i-1)*(length(y)):1+(i-1)*...
                (length(y))+(length(y)-1),1:length(y));
            dzydzij=[dzydzij_temp1(j_temp0,j)*(i==i_temp0)...
                dzydzij_temp1(j_temp0,j)*(i==i_temp0+1)...
                dzydzij_temp1(j_temp0+1,j)*(i==i_temp0)...
                dzydzij_temp1(j_temp0+1,j)*(i==i_temp0+1)];
            
            dzxydzij_temp1=Zxy(1+(j-1)*(length(x)):1+(j-1)*...
                (length(x))+(length(x)-1),1:length(x));
            dzxydzij=[dzxydzij_temp1(i_temp0,i)*(j==j_temp0)...
                dzxydzij_temp1(i_temp0+1,i)*(j==j_temp0)...
                dzxydzij_temp1(i_temp0,i)*(j==j_temp0+1)...
                dzxydzij_temp1(i_temp0+1,i)*(j==j_temp0+1)];%1x4
            
            
            if (i==i_temp0||i==i_temp0+1)&&(j==j_temp0||j==...
                    j_temp0+1)
                if i==i_temp0&&j==j_temp0
                    dgstdzij=dgstdzij_temp;
                elseif i==i_temp0+1&&j==j_temp0
                    dgstdzij=dgstdzi1j_temp;
                elseif i==i_temp0&&j==j_temp0+1
                    dgstdzij=dgstdzij1_temp;
                else
                    dgstdzij=dgstdzi1j1_temp;
                end
            else
                dgstdzij=0; 
            end

                dgstdzx=[dgstdzxij_temp,dgstdzxi1j_temp...
                    ,dgstdzxij1_temp,dgstdzxi1j1_temp];
                
                dgstdzy=[dgstdzyij_temp,dgstdzyi1j_temp...
                    ,dgstdzyij1_temp,dgstdzyi1j1_temp];
                
                dgstdzxy=[dgstdzxyij_temp,...
                    dgstdzxyi1j_temp,dgstdzxyij1_temp,...
                    dgstdzxyi1j1_temp]; %4x16
                
                dcstdzij=dgstdzij'+dzxdzij*dgstdzx'+...
                    dzydzij*dgstdzy'+dzxydzij*dgstdzxy';
                
                dfmdcst=xMyM(m,1:16); 
                
                sum1=dcstdzij*dfmdcst';

            p_n=xMyM(m,1:16)';
            c_n=c(1:16,n)';
            delzij_temp(m,1)=sum1*sign(c_n*p_n-zM(m));
            %             epsilon(m)=(c_n*p_n-zM(m)).^2;
            %             delzij_temp(m,1)=sum1*2*(c_n*p_n-zM(m));
        end
        %% Store gradient values for each spline node (i,j)
        delzij(i,j)=sum(delzij_temp);
    end
end
end