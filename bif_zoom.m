%% TXT_FILES

% 1 -> 3 -> 2 RE-ENTRANT BEHAVIOUR
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=3;b=1;
A_star_lower = 0.4;
A_star_upper = 0.46;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%starting time of loop
clock_start = datestr(now,'HH:MM:SS');

%fsolve tolerances
options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12);

k=1;n=4;theta_a=0.5;theta_b=0.5;

%pre-setting matrix size to speed up computations and setting inital row in matrix as zero
matrix_row = 0;steady_state_matrix=zeros(400000,13);
%energy values to scan through and step size
for A_star1 = A_star_lower:0.0005:A_star_upper
    A_star=round(A_star1,4);
    %display where the computation is up to in command window - good for long computations to see where you are up to
    fprintf('a = %.2f, b = %.2f, A*=%.4f at %s.\n',a,b,A_star,datestr(now,'HH:MM:SS'));
    %lambda
    l= @(A_star) 1./(1+exp(-(16*A_star-8)));
    %function to be used by fsolve
    fhandle=@(X)ode_sigmoid_model(X,a,b,k,A_star,n,theta_a,theta_b);
    
    %symbolic variables for protein levels
    syms x1 x2;
    %ODEs
    ode_eqns = [l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;...
        l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
    % variables for jacobian matrix
    ode_variables=[x1,x2];
    %calculating jacbian with respect to variables x1 & x2
    jac=jacobian(ode_eqns,ode_variables);
    
    for i=0:0.2:2%initial conditions on x axis
        for j=0:0.2:2%initial conditions on y axis
            
            %initial conditions pairing
            X0 = [i,j];
            %using fsolve
            [X,fval,exitflag,output] = fsolve(fhandle,X0,options);
            
            if (X(1) >= 0) && (X(2) >=0)
                %increasing matrix row by 1 in each loop
                matrix_row=matrix_row+1;
                
                %rounded calculated steady state values
                x1_ss=round(X(1),3);x2_ss=round(X(2),3);
                
                %subs. in steady state values to jacobian
                sub=subs(jac, [x1 x2], [X(1) X(2)]);
                %calc eigenvlaues
                eigen = eig(sub);
                %calculate the sign of each eigenvalue
                eigenvalue1=sign(eigen(1));eigenvalue2=sign(eigen(2));
                
                %testing if the steady state is stable or unstable
                if (eigenvalue1 < 0) && (eigenvalue2 < 0)
                    stability = 1;
                else
                    stability = -1;
                end
                
                %matrix of b, energy, lambda, ics, ss positions, the stability and subs2
                steady_state_matrix(matrix_row,:) = [a b A_star l(A_star) i j x1_ss x2_ss stability fval(1) fval(2) X(1) X(2)];
            else
                disp('negative steady state')
            end
        end
    end
end


txtFileName = sprintf('txt-files\\bif-n%d-a=%.0f-b=%.0f-zoom.txt',n,a*100,b*100);
fulltxtFileName=fullfile(txtFileName);
steady_state_matrix(~any(steady_state_matrix,2),:) = [];
fid = fopen(fulltxtFileName,'wt');
for ii = 1:size(steady_state_matrix,1)
    fprintf(fid,'%20.18f\t',steady_state_matrix(ii,:));
    fprintf(fid,'\n');
end

fprintf('b=%.2f a=%.2f zoomed bifurcation diagram. Start: %s. End: %s.\n',b,a,clock_start,datestr(now,'HH:MM:SS'));

% 1 -> 4 -> 3 RE-ENTRANT BEHAVIOUR
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1;b=0.5;
A_star_lower = 0.48;
A_star_upper = 0.54;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%starting time of loop
clock_start = datestr(now,'HH:MM:SS');

%fsolve tolerances
options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12);

k=1;n=4;theta_a=0.5;theta_b=0.5;

%pre-setting matrix size to speed up computations and setting inital row in matrix as zero
matrix_row = 0;steady_state_matrix=zeros(400000,13);
%energy values to scan through and step size
for A_star1 = A_star_lower:0.0005:A_star_upper
    A_star=round(A_star1,4);
    %display where the computation is up to in command window - good for long computations to see where you are up to
    fprintf('a = %.2f, b = %.2f, A*=%.4f at %s.\n',a,b,A_star,datestr(now,'HH:MM:SS'));
    %lambda
    l= @(A_star) 1./(1+exp(-(16*A_star-8)));
    %function to be used by fsolve
    fhandle=@(X)ode_sigmoid_model(X,a,b,k,A_star,n,theta_a,theta_b);
    
    %symbolic variables for protein levels
    syms x1 x2;
    %ODEs
    ode_eqns = [l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;...
        l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
    % variables for jacobian matrix
    ode_variables=[x1,x2];
    %calculating jacbian with respect to variables x1 & x2
    jac=jacobian(ode_eqns,ode_variables);
    
    for i=0:0.2:2%initial conditions on x axis
        for j=0:0.2:2%initial conditions on y axis
            
            %initial conditions pairing
            X0 = [i,j];
            %using fsolve
            [X,fval,exitflag,output] = fsolve(fhandle,X0,options);
            
            if (X(1) >= 0) && (X(2) >=0)
                %increasing matrix row by 1 in each loop
                matrix_row=matrix_row+1;
                
                %rounded calculated steady state values
                x1_ss=round(X(1),3);x2_ss=round(X(2),3);
                
                %subs. in steady state values to jacobian
                sub=subs(jac, [x1 x2], [X(1) X(2)]);
                %calc eigenvlaues
                eigen = eig(sub);
                %calculate the sign of each eigenvalue
                eigenvalue1=sign(eigen(1));eigenvalue2=sign(eigen(2));
                
                %testing if the steady state is stable or unstable
                if (eigenvalue1 < 0) && (eigenvalue2 < 0)
                    stability = 1;
                else
                    stability = -1;
                end
                
                %matrix of b, energy, lambda, ics, ss positions, the stability and subs2
                steady_state_matrix(matrix_row,:) = [a b A_star l(A_star) i j x1_ss x2_ss stability fval(1) fval(2) X(1) X(2)];
            else
                disp('negative steady state')
            end
        end
    end
end


txtFileName = sprintf('txt-files\\bif-n%d-a=%.0f-b=%.0f-zoom.txt',n,a*100,b*100);
fulltxtFileName=fullfile(txtFileName);
steady_state_matrix(~any(steady_state_matrix,2),:) = [];
fid = fopen(fulltxtFileName,'wt');
for ii = 1:size(steady_state_matrix,1)
    fprintf(fid,'%20.18f\t',steady_state_matrix(ii,:));
    fprintf(fid,'\n');
end

fprintf('b=%.2f a=%.2f zoomed bifurcation diagram. Start: %s. End: %s.\n',b,a,clock_start,datestr(now,'HH:MM:SS'));

%% FIGURES

% 1 -> 3 -> 2 RE-ENTRANT BEHAVIOUR
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=3;b=1;
A_star_lower = 0.4;
A_star_upper = 0.46;
n=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%font, fontsize and figure size
fn='Helvetica';wd=10;ht=9;fs_axis=11;%fs_labels=10;

txtFileName = sprintf('txt-files\\bif-n%d-a=%.0f-b=%.0f-zoom.txt',n,a*100,b*100);
steady_state_matrix = importdata(txtFileName);

%creating fig#ure
fig1=figure('Name','bifurcation_zoom');

%checking if the steady state when sub into ODEs is within a tolerance of plus-minus 10^-10
col_check=steady_state_matrix(:,10);ss_2=steady_state_matrix(abs(col_check) < 1e-8,:);
col_check2=ss_2(:,11);ss_5=ss_2(abs(col_check2) < 1e-8,:);

%extracting stable and unstable sub-matrices
col_stable1=ss_5(:,9);ss_1c=ss_5(col_stable1==1,:);ss_6=unique(ss_1c(:,[1 3 2 7 8]),'rows');
col_unstable1=ss_5(:,9);ss_1d=ss_5(col_unstable1==-1,:);ss_7=unique(ss_1d(:,[1 3 2 7 8]),'rows');

%plotting sub-figre
xlim([A_star_lower A_star_upper]);hold on;ylim([0 1.5]);grid on;box on;%axes limits & box around figure
plot(ss_6(:,2),ss_6(:,5),'b.','MarkerSize',4.5);%plotting stable steady states
plot(ss_7(:,2),ss_7(:,5),'r.','MarkerSize',4.5);%plotting unstable steady states
%hx=xlabel('$A^{*}$','interpreter','latex');hx.FontSize=fs_labels;hx.FontName=fn;%label on x-axis
%hy=ylabel('$x_{2}$-Position of Fixed Points','interpreter','latex');hy.FontSize=fs_labels;hy.FontName=fn;%label on y-axis
ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';ax.XTick = A_star_lower:0.01:A_star_upper;ax.YTick=0:0.25:1.5;%changing x and y axes properties
fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];%setting figure size
hold off;

%saving produced figure to output directory with specified name and file extenstion
epsFileName = sprintf('figures\\bif-n%d-a=%.0f-b=%.0f-zoom.eps',n,a*100,b*100);print(fig1,epsFileName,'-depsc');
tiffFileName = sprintf('figures\\bif-n%d-a=%.0f-b=%.0f-zoom.tiff',n,a*100,b*100);print(fig1,tiffFileName,'-dtiff');


% 1 -> 4 -> 3 RE-ENTRANT BEHAVIOUR
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1;b=0.5;
A_star_lower = 0.48;
A_star_upper = 0.54;
n=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%font, fontsize and figure size
fn='Helvetica';wd=10;ht=9;fs_axis=11;%fs_labels=10;

txtFileName = sprintf('txt-files\\bif-n%d-a=%.0f-b=%.0f-zoom.txt',n,a*100,b*100);
steady_state_matrix = importdata(txtFileName);

%creating fig#ure
fig1=figure('Name','bifurcation_zoom');

%checking if the steady state when sub into ODEs is within a tolerance of plus-minus 10^-10
col_check=steady_state_matrix(:,10);ss_2=steady_state_matrix(abs(col_check) < 1e-10,:);
col_check2=ss_2(:,11);ss_5=ss_2(abs(col_check2) < 1e-10,:);

%extracting stable and unstable sub-matrices
col_stable1=ss_5(:,9);ss_1c=ss_5(col_stable1==1,:);ss_6=unique(ss_1c(:,[1 3 2 7 8]),'rows');
col_unstable1=ss_5(:,9);ss_1d=ss_5(col_unstable1==-1,:);ss_7=unique(ss_1d(:,[1 3 2 7 8]),'rows');

%plotting figure
xlim([A_star_lower A_star_upper]);hold on;ylim([0 1]);grid on;box on;%axes limits & box around figure
plot(ss_6(:,2),ss_6(:,5),'b.','MarkerSize',4.5);%plotting stable steady states
plot(ss_7(:,2),ss_7(:,5),'r.','MarkerSize',4.5);%plotting unstable steady states
%hx=xlabel('$A^{*}$','interpreter','latex');hx.FontSize=fs_labels;hx.FontName=fn;%label on x-axis
%hy=ylabel('$x_{2}$-Position of Fixed Points','interpreter','latex');hy.FontSize=fs_labels;hy.FontName=fn;%label on y-axis
ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';ax.XTick = A_star_lower:0.01:A_star_upper;%changing x and y axes properties
fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];%setting figure size
hold off;

%saving produced figure to output directory with specified name and file extenstion
epsFileName = sprintf('figures\\bif-n%d-a=%.0f-b=%.0f-zoom.eps',n,a*100,b*100);print(fig1,epsFileName,'-depsc');
tiffFileName = sprintf('figures\\bif-n%d-a=%.0f-b=%.0f-zoom.tiff',n,a*100,b*100);print(fig1,tiffFileName,'-dtiff');