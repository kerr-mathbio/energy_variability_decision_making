%% ARRAYS FOR STABLE AND UNSTABLE POINTS WHEN CHANGING ATP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1;b=0.25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TXT FILES
k=1;n=4;theta_a=0.5;theta_b=0.5;

options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12);

m=0;
ss=zeros(5000,10);
stable_ss=zeros(2500,4);

for A_star = [0.4 0.5 0.6 0.7]
    l= @(A_star) 1./(1+exp(-(16*A_star-8)));%lambda([ATP])
    fprintf('Running array with a=%.2f, b=%.2f and A*=%.3f at %s.\n',a,b,A_star,datestr(now,'HH:MM:SS'));
    
    syms x1 x2;
    f1 = [l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;...
        l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
    v=[x1,x2];
    jac=jacobian(f1,v); %calculating jacbian with respect to variables x1 & x2
    
    fhandle=@(X)ode_sigmoid_model(X,a,b,k,A_star,n,theta_a,theta_b);
    
    for i=0:0.2:2
        for j=0:0.2:2
            
            X0 = [i,j];
            [X,fval,exitflag,output] = fsolve(fhandle,X0,options);
            
            if (X(1) >=0) && (X(2) >=0)
                m=m+1;
                x1_ss=round(X(1),2);x2_ss=round(X(2),2);
                sub=subs(jac, [x1 x2], [X(1) X(2)]);
                eigen = eig(sub); %calc eigenvlaues of matrix 'sub'
                eigenvalue_1_2=sign(eigen(1));eigenvalue_2_2=sign(eigen(2));
                
                if (eigenvalue_1_2 < 0) && (eigenvalue_2_2 < 0)
                    stability = 1;
                else
                    stability = -1;
                end
                %matrix of value of a, b, energy, ics, ss positions for x1, x2, the stability of the ss point and steady states evaluated in ODEs
                ss(m,:) = [a b A_star i j x1_ss x2_ss stability fval(1) fval(2)];
            else
                disp('negative ss');
            end
        end
    end
end

txtFileName = sprintf('txt-files\\array-atp-a=%.0f-b=%.0f.txt',a*100,b*100);
fulltxtFileName=fullfile(txtFileName);
ss(~any(ss,2),:) = [];
fid = fopen(fulltxtFileName,'wt');
for ii = 1:size(ss,1)
    fprintf(fid,'%20.18f\t',ss(ii,:));
    fprintf(fid,'\n');
end

%% FIGURES

fn='Helvetica';wd=16;ht=14;fs_labels=12;fs_axis=12;%fs_title=12;

k=1;n=4;theta_a=0.5;theta_b=0.5;

txtFileName = sprintf('txt-files\\array-atp-a=%.0f-b=%.0f.txt',a*100,b*100);
ss = importdata(txtFileName);

col_check=ss(:,9);ss_b=ss(col_check < 1e-8,:);
col_check2=ss_b(:,10);ss_c=ss_b(col_check2 < 1e-8,:);
col_check3=ss_c(:,9);ss_d=ss_c(col_check3 > -1e-8,:);
col_check4=ss_d(:,10);ss_1=ss_d(col_check4 > -1e-8,:);

array_fig=figure;
set(gcf,'Units','centimeters','Position',[0 0 wd ht],'PaperUnits','centimeters','PaperSize',[wd ht]);
p=0;

for A_star = [0.4 0.5 0.6 0.7]
    p=p+1;
    
    col_a1=ss_1(:,1);
    ss_1b = ss_1(col_a1 == a,:);
    col_b1=ss_1b(:,2);
    ss_1c = ss_1b(col_b1 == b,:);%new matrix B_2 is a submatrix of B with value a in the first column
    col_y1 = ss_1c(:,3);%selecting a column in B_2 matrix
    ss_1d = ss_1c(col_y1 == A_star,:);%new matrix B_3 is a submatrix of B_2 with value b in the second column
    
    stab_col1=ss_1d(:,8);stable_ss1 = ss_1d(stab_col1==1,:);unique_stable_ss1 = unique(stable_ss1(:,[1 2 3 6 7]),'rows');
    unstab_col1=ss_1d(:,8);unstable_ss1 = ss_1d(stab_col1==-1,:);unique_unstable_ss1 = unique(unstable_ss1(:,[1 2 3 6 7]),'rows');
    
    
    subplot(2,2,p);
    plot(unique_unstable_ss1(:,4),unique_unstable_ss1(:,5),'ro','MarkerSize',3,'MarkerFaceColor','r');hold on;
    plot(unique_stable_ss1(:,4),unique_stable_ss1(:,5),'bo','MarkerSize',3,'MarkerFaceColor','b');hold off;%,'MarkerSize',2);
    grid on;box on; %axes limits & box around figure
    ax = gca;ax.XTick = 0:1:4;ax.YTick = 0:1:4;ax.YLim = [0 4];ax.XLim = [0 4];
    ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';
% %     hx=xlabel('$x_1$');hx.Interpreter='latex';hx.FontSize=fs_labels;hx.FontName=fn;%x-axis
% %     hy=ylabel('$x_2$');hy.Interpreter='latex';hy.FontSize=fs_labels;hy.FontName=fn;%y-axis
% %     sub_tit=title(sprintf('A* = %.4f',A_star)); %changing title within for loop for value of ATP used in calculations
% %     sub_tit.FontName=fn;sub_tit.FontWeight='normal';sub_tit.FontSize=12;%fs_axis;
    %   lgd=legend('Steady State Points','Stable Steady States');%,lgd.Position=[left bottom (thse two do dist from lower left corner to lower left corner of legend) width height (legend dimensions)]; best adding units here too lgd.Units='centimeters';
    %   lgd.FontSize=legend_fs;lgd.FontName=fn;lgd.Interpreter='LaTex';lgd.Location='northeast'; %legend properties
    
   
    % %                 ax = gca;ax.XTick = 0:1:6;ax.YTick = 0:2:6;ax.YLim = [0 6];ax.XLim = [0 6];
end

epsFileName = sprintf('figures\\array-atp-a=%.0f-b=%.0f.eps',a*100,b*100);fullFileName =fullfile(epsFileName);print(array_fig,fullFileName,'-depsc');
tiffFileName = sprintf('figures\\array_atp-a=%.0f-b=%.0f.tiff',a*100,b*100);fullFileName=fullfile(tiffFileName);print(array_fig,fullFileName,'-dtiff');