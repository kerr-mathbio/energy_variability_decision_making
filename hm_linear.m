%% HEATMAPS FOR LINEAR LAMBDA

%% TXT_FILES
clear;
%set a values to produce heatmaps for
for a=[0,0.5,1,1.5,2,3]
    
    %fsolve tolerances
    options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12);
    
    %fixed parameter values for ODEs
    k=1;n=4;theta_a=0.5;theta_b=0.5;
    
    %pre-setting matrix size to speed up computations
    ss_matrix=zeros(1200000,11);
    %setting the initial row in matrix 'ss_matrix' as zero
    matrix_row=0;
    
    %starting time of loop
    clock_start = datestr(now,'HH:MM:SS');
    
    %set b values
    for b=0:0.25:3
        %energy values to scan through
        for A_star = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]
            %lambda
            l= @(A_star) A_star;
            %display where the computation is up to in command window - good for long computations to see where you are up to
            fprintf('Running heatmap with a=%.2f, b=%.2f and A_star=%.2f at %s.\n',a,b,A_star,datestr(now,'HH:MM:SS'));
            %fsolve function
            fhandle=@(X)ode_linear_model(X,a,b,k,A_star,n,theta_a,theta_b);
            %symbolic variables for protein levels
            syms x1 x2;
            %ODEs
            f1 = [l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;...
                l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
            %variables for jacobian matrix
            v=[x1,x2];
            %calculating jacbian with respect to variables x1 & x2
            jac=jacobian(f1,v);
            
            %initial conditions on x axis
            for i=0:0.2:2
                %initial conditions on y axis
                for j=0:0.2:2
                    %initial conditions pairing
                    X0 = [i,j];
                    %using fsolve
                    [X,fval,exitflag,output] = fsolve(fhandle,X0,options);
                    
                    if (X(1) >=0) && (X(2) >=0)
                        %moving row in matrix down by 1 in each loop
                        matrix_row=matrix_row+1;
                        
                        %rounded calculated steady state values
                        x1_ss=round(X(1),3);x2_ss=round(X(2),3);
                        
                        %subs. in steady state values to jacobian
                        sub=subs(jac, [x1 x2], [X(1) X(2)]); %subs. in ss values from original ics
                        %calc eigenvlaues
                        eigen = eig(sub); %calc eigenvlaues of matrix 'sub'
                        %calculate the sign of each eigenvalue
                        eigenvalue_1=sign(eigen(1));eigenvalue_2=sign(eigen(2));
                        
                        %testing if the steady state is stable or unstable
                        if (eigenvalue_1 < 0) && (eigenvalue_2 < 0)
                            stability = 1;
                        else
                            stability = -1;
                        end
                        
                        %rounding lambda
                        lambda = round(l(A_star),3);
                        %matrix of a, b, energy, ics, ss positions, the stability, lambda and subs2
                        ss_matrix(matrix_row,:) = [a b A_star i j x1_ss x2_ss stability lambda fval(1) fval(2)];
                    else
                        disp('negative steady state');
                    end
                end
            end
        end
    end
    
    txtFileName = sprintf('txt-files\\hm-linear-n%d-a=%.0f.txt',n,a*100);
    fulltxtFileName=fullfile(txtFileName);
    ss_matrix(~any(ss_matrix,2),:) = [];
    fid = fopen(fulltxtFileName,'wt');
    for ii = 1:size(ss_matrix,1)
        fprintf(fid,'%20.18f\t',ss_matrix(ii,:));
        fprintf(fid,'\n');
    end
    
    fprintf('\n \n Linear: a=%.2f, b = %.2f, A*=%.2f. Start: %s. End: %s\n \n',a,b,A_star,clock_start,datestr(now,'HH:MM:SS'));
    
    clear;
    
    pause(30);
    
end


%% FIGURES

%universal values
wd=10;ht=9;
fn='Helvetica';
fs_labels=12;
fs_axis=11;

%order of energy concetrations for y-axis -- see use further down
NO = {'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'};

%set activation values to produce heatmaps for
for a=[0,0.5,1,1.5,2,3]
    
    n=4;
    
    txtFileName = sprintf('txt-files\\hm-linear-n%d-a=%.0f.txt',n,a*100);
    ss_matrix = importdata(txtFileName);
    
    %checking if the steady state when sub into ODEs is within a tolerance of plus-minus 10^-10
    col_check=ss_matrix(:,10);ss_matrix2=ss_matrix(abs(col_check) < 1e-8,:);
    col_check2=ss_matrix2(:,11);ss_1=ss_matrix2(abs(col_check2) < 1e-8,:);
    
    %selecting a column in ss_1 matrix
    stab_col3=ss_1(:,8);
    %new sub-matrix ss_8 is a submatrix of ss_1 with stable steady states
    ss_8=ss_1(stab_col3 == 1,:);
    %extracting unique stable steady states
    ss_9 = unique(ss_8(:,[1 2 3 6 7 9]),'rows');
    %selecting a column in ss_9 matrix
    col_a4 = ss_9(:,1);
    %new sub-matrix stable_ss3 is a submatrix of ss_9 with a=a
    stable_ss3 = ss_9(col_a4 == a,:);
    
    %table of a, b, energy, stable steady state positions and lambda value
    T2 = array2table(stable_ss3,'VariableNames',{'a','b','energy','stable_ss_position_x1','stable_ss_position_x2','lambda'});
    
    %creating figure
    fig_heatmap = figure('Name','Heatmap');
    %plotting heatmap from table T2' with b on the x-axis and energy values on the y-axis
    h = heatmap(T2,'b','energy','Title','');
    %colour scheme to use -- could use spring or summer or parula
    h.Colormap = cool;
    %data colour and label if missing
    h.MissingDataColor = [0.8 0.8 0.8];h.MissingDataLabel = 'No data';
    %method to use for displaying data
    h.ColorMethod = 'count';
    %color bar visible or not in figures
    h.ColorbarVisible = 'off';
    %fontname, colour and fontsize
    h.FontColor = 'black';h.FontName = fn;h.FontSize = fs_labels;
    %colour limits to keep consistency between figures
    h.ColorLimits=[1 4];
    
    %re-arranging y-axis energy values
    h.SourceTable.energy = categorical(h.SourceTable.energy);
    neworder = NO;
    h.SourceTable.energy = reordercats(h.SourceTable.energy,neworder);
    
    %removing y-axis label
    h.YLabel=' ';
    %removing x-axis label
    h.XLabel=' ';
    %axis fontname and fontsize
    ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;
    %figure size
    fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];
    
    %saving produced figure to output directory with specified name and file extenstion
    epsFileName = sprintf('figures\\hm-linear-n%d-a=%.0f.eps',n,a*100);fullFileName = fullfile(epsFileName);print(fig_heatmap,fullFileName,'-depsc');
    tiffFileName = sprintf('figures\\hm-linear-n%d-a=%.0f.tiff',n,a*100);fullFileName2=fullfile(tiffFileName);print(fig_heatmap,fullFileName2,'-dtiff');
end
