%% Approximation Test
clear
clc
%% Data Acquisition
A=[652 725 102 731 506 79 223 438 767 772 652 725 102 731 506 79 223 438 767 772 652 725 102 731 506 79 223 438 767 772 652 725 102 731 506 79 223 438 767 772 652 725 102 731 506 79 223 438 767 772 652 725 102 731 506 79 223 438 767 772];%Killeen_100x100
% Define input data size with N where N = n means (n+1) x (n+1) dataset
n = 9;
C = A + n;
%Load Dataset
load('Killeen_901x901.mat');
%Choose Subset Range
for counter=1
    x1=A(counter);
    x2=C(counter);
    y1=x1;
    y2=x2;
    zMtemp=double(J(x1:x2,y1:y2,:));    
    %% Pre-processing
    xM_given=(1:length(x1:x2));
    yM_given=(1:length(y1:y2));
    [Xtemp,Ytemp]=meshgrid(xM_given,yM_given);
    xM=reshape(Xtemp',1,numel(Xtemp));
    yM=reshape(Ytemp',1,numel(Ytemp));
    zM=reshape(zMtemp',1,numel(zMtemp));
    % Define ratio of given data to node size where ratio = r means r x r
    % nodes are selected.
    ratio = 0.5;
    x=linspace(1,n+1,(n+1)*ratio);
    y=x;
    % Calculate constant xMyM matrix
    xMyM=[];
    I=find(x==max(x));
    J1=find(y==max(y));
    for m=1:length(xM)
        [i_temp0,j_temp0]=Mapping_Phase_1(x,y,xM(m),yM(m));
        [~,xMyM_temp]=bicubic_grid((xM(m)-x(i_temp0))/(x(i_temp0+1)-x(i_temp0))...
            ,(yM(m)-y(j_temp0))/(y(j_temp0+1)-y(j_temp0)));
        xMyM=cat(1,xMyM,xMyM_temp);
    end
    %% Initialization
    tic
    z0=zeros();
    X = [ones(size(xM')) xM' yM' xM'.*yM'];
    b = regress(zM',X);
    for i=1:length(x)
        for j=1:length(y)
            z0(i,j)=b(1) + b(2).*x(i) + b(3).*y(j) + b(4).*x(i).*y(j);
        end
    end
    %% Backtracking Strategy
    alpha=0.25;
    beta=0.5;
    s=2;
    epsilon=10^-3;
    z=z0;
    delzij=Gradient_finder(x,y,z0,xM,yM,zM,xMyM);
    c=Bicubic_interpolation(x,y,z0);
    fun_val=fiterr(xM,yM,zM,xMyM,x,y,I,J1,c);
    Iteration_number=zeros();
    Apprx_error=zeros();
    t=1;
    iter=0;
    while iter<100
        iter=iter+1;
        t=s;
        c=Bicubic_interpolation(x,y,z-(t*delzij/norm(delzij,'fro')));
        f=fiterr(xM,yM,zM,xMyM,x,y,I,J1,c);
        while ((fun_val-f)<alpha*t&&t>epsilon)
            t=beta*t;
            c=Bicubic_interpolation(x,y,z-(t*delzij/norm(delzij,'fro')));
            f=fiterr(xM,yM,zM,xMyM,x,y,I,J1,c);
        end
        Iteration_number(iter,1)=iter;
        z=z-(t*delzij/norm(delzij,'fro'));
        delzij=Gradient_finder(x,y,z,xM,yM,zM,xMyM);
        c=Bicubic_interpolation(x,y,z);
        fun_val=fiterr(xM,yM,zM,xMyM,x,y,I,J1,c);
        Apprx_error(iter,1)=fun_val/(size(x,2)*size(y,2));
        if iter>1
            if (1-(Apprx_error(iter,1))/Apprx_error(iter-1,1))<0.001 &&...
                    (1-(Apprx_error(iter,1))/Apprx_error(iter-1,1))>0
                break
            end
        end
        fprintf('Iteration no = %d Norm of Gradient = %2.6f Approximation Error = %2.6f \n',...
            iter,norm(delzij,'fro'),fun_val/(size(xM,2)*size(yM,2)));
    end
    e1=toc;
    disp('Optimization Terminated')
    
    z_final=z;
    iteration=iter;
    %% Plot Parameters
    
    % 4 plots will be generated as follows:
    % 1. Error Minimization Trend
    % 2. Approximated 3D Spline Surface
    % 3. Comparison of Given, Interpolated and Approximated Spline Surfaces
    % 4. Comparison of Given, Interpolated and Approximated surface contour
    
    
    % 1. Error Minimization Trend plot
    
    % Here the error minimization trend with each iteration is displayed
    figure('Name','Iteration Number by Approximation Error')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    plot(Iteration_number,Apprx_error,'r-o')
%     print(strcat('ItervsError_',num2str(counter)),'-djpeg')
    
    % 2. Approximated 3D Spline Surface plot
    
    % Here the interpolation is performed with the optimized nodes
    [xplot,yplot,zplot]=drawspline2D(x,y,z_final,10,10);
    figure('Name','Approximated 2-D Spline Surface')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    surf(xplot,yplot,zplot');
    hold on
    scatter3(xM,yM,zM,'+')
    hold on
    [X,Y]=meshgrid(x,y);
    scatter3(Y(:),X(:),z_final(:),100,'filled')
    hold off
    legend('Spline Surface',['Given Data Points [' num2str(size(zMtemp)) ']'],...
        ['Optimized Nodes [' num2str(size(z_final)) ']']')
    set(legend,'FontSize',12,'Location','best')
%     print(strcat('ApproximatedSpline_',num2str(counter)),'-djpeg')
%     close all
    
    %3. Comparison of Given, Interpolated and Approximated Spline Surfaces
    
    % Here, the interpolation is performed with each data point as a node
    [xi,yi,zi]=drawspline2D(xM_given,yM_given,zMtemp,20,20);
    
    
    figure('Name','Given, Interpolated and Approximated 2-D Spline Surface')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    % Equal Axes for each subplot
    
    % I couldn't do it easily but I am sure there is an easier way
    
    
    m1=max(max((zi)));
    m2=max(max(zplot));
    if m1>m2
        m_com=m1;
    else
        m_com=m2;
    end
    n1=min(min((zi)));
    n2=min(min(zplot));
    if n1<n2
        n_com=n1;
    else
        n_com=n2;
    end
    subplot (1,3,1),surf(xM_given, yM_given, zMtemp)
    shading interp
    title('Given','FontSize',12)
    axis([-inf inf -inf inf floor(n_com) ceil(m_com)])
    subplot (1,3,2),surf(xi,yi,zi)
    shading interp
    title('Interpolation','FontSize',12)
    axis([-inf inf -inf inf floor(n_com) ceil(m_com)])
    subplot (1,3,3),surf(xplot,yplot,zplot')
    shading interp
    title('Approximation','FontSize',12)
    axis([-inf inf -inf inf floor(n_com) ceil(m_com)])
%     print(strcat('Given, Interpolated and Approximated 2-D Spline Surfaces'...
%         ,num2str(counter)),'-djpeg')
%     close all
    
    
    % 4. Comparison of Given, Interpolated and Approximated surface contour
    
    figure('Name','Given, Interpolated and Approximated Surface Contours')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
    subplot (1,3,1),contour(xM_given, yM_given, zMtemp)
    title('Given','FontSize',12)
    subplot (1,3,2),contour(xi,yi,zi)
    title('Interpolation','FontSize',12)
    subplot (1,3,3),contour(xplot,yplot,zplot')
    
    title('Approximation','FontSize',12)
%     print(strcat('Given, Interpolated and Approximated Surface Contours',...
%         num2str(counter)),'-djpeg')
%     close all
end
