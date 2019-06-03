%% BIFURCATION DIAGRAMS FOR a, b AND k

%% TXT-FILES

%BIFURCATION DIAGRAMS WHEN CHANGING b
clear;
%b values to scan through
for b=[0.25,0.5,0.75,1.5]
    
    %fixed parameter values for ODEs in these bifurcation diagrams
    a=1;k=1;n=4;theta_a=0.5;theta_b=0.5;
    
    %fsolve tolerances
    options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12);
    
    %pre-setting matrix size to speed up computations and setting inital row in matrix as zero
    matrix_row1 = 0;ss_1a=zeros(55000,10);
    
    %starting time of loop
    clock_start = datestr(now,'HH:MM:SS');
    
    %energy values to scan through and step size
    for A_star = 0:0.005:1.0%0:0.05:1.0%0:0.1:1%
        %display where the computation is up to in command window - good for long computations to see where you are up to
        fprintf('b=%.2f bifurcation diagram, A_star=%.4f at %s.\n',b,A_star,datestr(now,'HH:MM:SS'));
        %lambda
        l1= @(A_star) 1./(1+exp(-(16*A_star-8)));%lambda
        %function to be used by fsolve
        fhandle=@(X)ode_sigmoid_model(X,a,b,k,A_star,n,theta_a,theta_b);
        
        %symbolic variables for protein levels
        syms x1 x2;
        %ODEs
        f1 = [l1(A_star)*a*x1^n./(theta_a^n+x1^n)+l1(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;...
            l1(A_star)*a*x2^n./(theta_a^n+x2^n)+l1(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
        % variables for jacobian matrix
        v1=[x1,x2];
        %calculating jacbian with respect to variables x1 & x2
        jac1=jacobian(f1,v1);
        
        for i=0:0.2:2%initial conditions on x axis
            for j=0:0.2:2%initial conditions on y axis
                
                %initial conditions pairing
                X0 = [i,j];
                %using fsolve
                [X,fval,exitflag,output] = fsolve(fhandle,X0,options);
                
                if (X(1) >= 0) && (X(2) >=0)
                    %increasing matrix row by 1 in each loop
                    matrix_row1=matrix_row1+1;
                    
                    %rounded calculated steady state values
                    x1_ss_1=round(X(1),3);x2_ss_1=round(X(2),3);
                    
                    %subs. in steady state values to jacobian
                    sub1=subs(jac1, [x1 x2], [X(1) X(2)]);
                    %calc eigenvlaues
                    eigen1 = eig(sub1);
                    %calculate the sign of each eigenvalue
                    eigenvalue_1_a=sign(eigen1(1));eigenvalue_2_a=sign(eigen1(2));
                    
                    %testing if the steady state is stable or unstable
                    if (eigenvalue_1_a < 0) && (eigenvalue_2_a < 0)
                        stability = 1;
                    else
                        stability = -1;
                    end
                    
                    %matrix of b, energy, lambda, ics, ss positions, the stability and subs2
                    ss_1a(matrix_row1,:) = [b A_star l1(A_star) i j x1_ss_1 x2_ss_1 stability fval(1) fval(2)];
                else
                    disp(' ')
                end
            end
        end
    end
    
    txtFileName = sprintf('txt-files\\bif-n%d-b=%.0f.txt',n,b*100);
    fulltxtFileName=fullfile(txtFileName);
    ss_1a(~any(ss_1a,2),:) = [];
    fid = fopen(fulltxtFileName,'wt');
    for ii = 1:size(ss_1a,1)
        fprintf(fid,'%20.18f\t',ss_1a(ii,:));
        fprintf(fid,'\n');
    end
    
    fprintf('b=%.2f bifurcation diagram. Start: %s. End: %s.\n',b,clock_start,datestr(now,'HH:MM:SS'));
    
    clear;
    
    pause(30);
    
end


%BIFURCATION DIAGRAMS WHEN CHANGING a
clear;
for a=0%[0,1,1.25,1.75,3]
    
    %fixed parameter values for ODEs in these bifurcation diagrams
    b=1;k=1;n=4;theta_a=0.5;theta_b=0.5;
    
    %fsolve tolerances
    options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12);
    
    matrix_row2=0;ss_2a=zeros(55000,10);
    
    clock_start = datestr(now,'HH:MM:SS');
    
    for A_star = 0:0.005:1.0
        fprintf('a=%.2f bifurcation diagram, A_star=%.4f at %s.\n',a,A_star,datestr(now,'HH:MM:SS'));
        l2= @(A_star) 1./(1+exp(-(16*A_star-8)));
        
        fhandle=@(X)ode_sigmoid_model(X,a,b,k,A_star,n,theta_a,theta_b);
        
        syms x1 x2;
        f2 = [l2(A_star)*a*x1^n./(theta_a^n+x1^n)+l2(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;...
            l2(A_star)*a*x2^n./(theta_a^n+x2^n)+l2(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
        v2=[x1,x2];
        jac2=jacobian(f2,v2);
        
        for i=0:0.2:2
            for j=0:0.2:2
                
                X0 = [i,j];
                [X,fval,exitflag,output] = fsolve(fhandle,X0,options);
                
                if (X(1) >= 0) && (X(2) >=0)
                    matrix_row2=matrix_row2+1;
                    
                    x1_ss_2=round(X(1),3);x2_ss_2=round(X(2),3);
                    
                    sub2=subs(jac2, [x1 x2], [X(1) X(2)]);
                    eigen2 = eig(sub2);
                    eigenvalue_1_2a=sign(eigen2(1));eigenvalue_2_2a=sign(eigen2(2));
                    
                    if (eigenvalue_1_2a < 0) && (eigenvalue_2_2a < 0)
                        stability = 1;
                    else
                        stability = -1;
                    end
                    
                    ss_2a(matrix_row2,:) = [a A_star l2(A_star) i j x1_ss_2 x2_ss_2 stability fval(1) fval(2)];
                else
                    disp(' ')
                end
            end
        end
    end
    
    txtFileName = sprintf('txt-files\\bif-n%d-a=%.0f.txt',n,a*100);
    fulltxtFileName=fullfile(txtFileName);
    ss_2a(~any(ss_2a,2),:) = [];
    fid = fopen(fulltxtFileName,'wt');
    for ii = 1:size(ss_2a,1)
        fprintf(fid,'%20.18f\t',ss_2a(ii,:));
        fprintf(fid,'\n');
    end
    
    fprintf('a=%.2f bifurcation diagram. Start: %s. End: %s.\n',a,clock_start,datestr(now,'HH:MM:SS'));
    
    clear;
    
    pause(30);
    
end

%BIFURCATION DIAGRAMS WHEN CHANGING k
clear;
for k=[0.5,1.25,1.5,3]
    
    %fixed parameter values for ODEs in these bifurcation diagrams
    a=1;b=1;n=4;theta_a=0.5;theta_b=0.5;
    
    %fsolve tolerances
    options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12);
    
    matrix_row3=0;ss_3a=zeros(55000,10);
    
    clock_start = datestr(now,'HH:MM:SS');
    
    for A_star1 = 0:0.005:1.0
        A_star=round(A_star1,3);
        fprintf('k=%.2f bifurcation diagram, A_star=%.4f at %s.\n',k,A_star,datestr(now,'HH:MM:SS'));
        l3= @(A_star) 1./(1+exp(-(16*A_star-8)));
        
        fhandle=@(X)ode_sigmoid_model(X,a,b,k,A_star,n,theta_a,theta_b);
        
        syms x1 x2;
        f3 = [l3(A_star)*a*x1^n./(theta_a^n+x1^n)+l3(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;...
            l3(A_star)*a*x2^n./(theta_a^n+x2^n)+l3(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
        v3=[x1,x2];
        jac3=jacobian(f3,v3);
        
        for i=0:0.2:2
            for j=0:0.2:2
                
                X0 = [i,j];
                [X,fval,exitflag,output] = fsolve(fhandle,X0,options);
                
                if (X(1) >= 0) && (X(2) >=0)
                    matrix_row3=matrix_row3+1;
                    x1_ss_3=round(X(1),3);x2_ss_3=round(X(2),3);
                    
                    sub3=subs(jac3, [x1 x2], [X(1) X(2)]);
                    eigen3 = eig(sub3);
                    eigenvalue_1_3a=sign(eigen3(1));eigenvalue_2_3a=sign(eigen3(2));
                    
                    if (eigenvalue_1_3a < 0) && (eigenvalue_2_3a < 0)
                        stability = 1;
                    else
                        stability = -1;
                    end
                    ss_3a(matrix_row3,:) = [k A_star l3(A_star) i j x1_ss_3 x2_ss_3 stability fval(1) fval(2)];
                else
                    disp(' ')
                end
            end
        end
    end
    
    txtFileName = sprintf('txt-files\\bif-n%d-k=%.0f.txt',n,k*100);
    fulltxtFileName=fullfile(txtFileName);
    ss_3a(~any(ss_3a,2),:) = [];
    fid = fopen(fulltxtFileName,'wt');
    for ii = 1:size(ss_3a,1)
        fprintf(fid,'%20.18f\t',ss_3a(ii,:));
        fprintf(fid,'\n');
    end
    
    fprintf('k=%.2f bifurcation diagram. Start: %s. End: %s.\n',k,clock_start,datestr(now,'HH:MM:SS'));
    
    clear;
    
    pause(30);
    
end


%% FIGURES

%universal figure properties
fn='Helvetica';wd=10;ht=9;fs_axis=12;%fs_labels=11;

% BIFURCATION DIAGRAMS WHEN CHANGING b   
%b values to scan through
for b=[0.25,0.5,0.75,1.5]
    
    txtFileName = sprintf('txt-files\\bif-n%d-b=%.0f.txt',n,b*100);
    ss_1a = importdata(txtFileName);

    %creating bifurcation figure
    fig1=figure('Name','Bifurcation_b');
    
    %checking if the steady state when sub into ODEs is within a tolerance of plus-minus 10^-10
    col_check=ss_1a(:,9);ss_1a2=ss_1a(abs(col_check) < 1e-8,:);
    col_check2=ss_1a2(:,10);ss_1a5=ss_1a2(abs(col_check2) < 1e-8,:);
   
    col_b=ss_1a5(:,1);%selecting a column in ss_1a5 matrix
    ss_1b=ss_1a5(col_b==b,:);%new sub-matrix ss_1b is a submatrix of ss_1a5 with value b=b in the first column
    
    %extracting stable and unstable sub-matrices
    col_stable1=ss_1b(:,8);ss_1c=ss_1b(col_stable1==1,:);ss_1e=unique(ss_1c(:,[1 2 6 7]),'rows');
    col_unstable1=ss_1b(:,8);ss_1d=ss_1b(col_unstable1==-1,:);ss_1f=unique(ss_1d(:,[1 2 6 7]),'rows');
    
    %plotting sub-figre
    xlim([0 1]);hold on;ylim([0 3]);grid on;box on;%axes limits & box around figure
    plot(ss_1e(:,2),ss_1e(:,4),'b.','MarkerSize',3);%plotting stable steady states
    plot(ss_1f(:,2),ss_1f(:,4),'r.','MarkerSize',3);%plotting unstable steady states
    %hx=xlabel('$A^{*}$','interpreter','latex');hx.FontSize=fs_labels;hx.FontName=fn;%label on x-axis
    %hy=ylabel('$x_{2}$-Position of Fixed Points','interpreter','latex');hy.FontSize=fs_labels;hy.FontName=fn;%label on y-axis
    ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';ax.XTick = 0:0.1:1;%changing x and y axes properties
    fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];%setting figure size
    hold off;
    
    %saving produced figure to output directory with specified name and file extenstion   
    epsFileName = sprintf('figures\\bif-n%d-b=%.0f.eps',n,b*100);fullFileName = fullfile(epsFileName);print(fig1,fullFileName,'-depsc');
    tiffFileName = sprintf('figures\\bif-n%d-b=%.0f.tiff',n,b*100);fullFileName2=fullfile(tiffFileName);print(fig1,fullFileName2,'-dtiff');
    
end

%clear some of the information stored by matlab - valuable if full code file is executed.
param={'b','a','k','y','i','j','x1','x2','X','X0','stability'};clear(param{:});


% BIFURCATION DIAGRAMS WHEN CHANGING a
for a=[0,1,1.25,1.75,3]
    
    txtFileName = sprintf('txt-files\\bif-n%d-a=%.0f.txt',n,a*100);
    ss_2a = importdata(txtFileName);

    fig2=figure('Name','Bifurcation_a');
    
    col_check=ss_2a(:,9);ss_2a2=ss_2a(abs(col_check) < 1e-8,:);
    col_check2=ss_2a2(:,10);ss_2a5=ss_2a2(abs(col_check2) < 1e-8,:);
    
    col_a=ss_2a5(:,1);ss_2b=ss_2a5(col_a==a,:);
    
    col_stable2=ss_2b(:,8);ss_2c=ss_2b(col_stable2==1,:);ss_2e=unique(ss_2c(:,[1 2 6 7]),'rows');
    col_unstable2=ss_2b(:,8);ss_2d=ss_2b(col_unstable2==-1,:);ss_2f=unique(ss_2d(:,[1 2 6 7]),'rows');
    
    xlim([0 1]);hold on;ylim([0 4]);grid on;box on;
    plot(ss_2e(:,2),ss_2e(:,4),'b.','MarkerSize',3);
    plot(ss_2f(:,2),ss_2f(:,4),'r.','MarkerSize',3);
    %hx=xlabel('$A^{*}$','interpreter','latex');hx.FontSize=fs_labels;hx.FontName=fn;
    %hy=ylabel('$x_{2}$-Position of Fixed Points','interpreter','latex');hy.FontSize=fs_labels;hy.FontName=fn;
    ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';ax.XTick = 0:0.1:1.0;
    fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];
    hold off;

    epsFileName = sprintf('figures\\bif-n%d-a=%.0f.eps',n,a*100);fullFileName = fullfile(epsFileName);print(fig2,fullFileName,'-depsc');
    tiffFileName = sprintf('figures\\bif-n%d-a=%.0f.tiff',n,a*100);fullFileName2=fullfile(tiffFileName);print(fig2,fullFileName2,'-dtiff');
    
    
end

%clear some of the information stored by matlab - valuable if full code file is executed.
param={'b','a','k','y','i','j','x1','x2','X','X0','stability'};clear(param{:});

% BIFURCATION DIAGRAMS WHEN CHANGING k
for k=[0.5,1.25,1.5,3]
    
    txtFileName = sprintf('txt-files\\bif-n%d-k=%.0f.txt',n,k*100);
    ss_3a = importdata(txtFileName);

    fig3=figure('Name','Degradation BD');
    
    col_check=ss_3a(:,9);ss_3a2=ss_3a(abs(col_check) < 1e-8,:);
    col_check2=ss_3a2(:,10);ss_3a5=ss_3a2(abs(col_check2) < 1e-8,:);
    
    col_k=ss_3a5(:,1);ss_3b=ss_3a5(col_k==k,:);
    col_stable3=ss_3b(:,8);ss_3c=ss_3b(col_stable3==1,:);ss_3e=unique(ss_3c(:,[1 2 6 7]),'rows');
    col_unstable3=ss_3b(:,8);ss_3d=ss_3b(col_unstable3==-1,:);ss_3f=unique(ss_3d(:,[1 2 6 7]),'rows');
    
    xlim([0 1]);hold on;ylim([0 4]);grid on;box on;
    plot(ss_3e(:,2),ss_3e(:,4),'b.','MarkerSize',3);
    plot(ss_3f(:,2),ss_3f(:,4),'r.','MarkerSize',3);
    %hx=xlabel('$A^{*}$','interpreter','latex');hx.FontSize=fs_labels;hx.FontName=fn;
    %hy=ylabel('$x_{2}$-Position of Fixed Points','interpreter','latex');hy.FontSize=fs_labels;hy.FontName=fn;
    ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';ax.XTick = 0:0.1:1.0;
    fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];
    hold off;

    epsFileName = sprintf('figures\\bif-n%d-k=%.0f.eps',n,k*100);fullFileName = fullfile(epsFileName);print(fig3,fullFileName,'-depsc');
    tiffFileName = sprintf('figures\\bif-n%d-k=%.0f.tiff',n,k*100);fullFileName2=fullfile(tiffFileName);print(fig3,fullFileName2,'-dtiff');
end
