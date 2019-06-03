%% ARRAYS FOR STABLE AND UNSTABLE POINTS WHEN CHANGING b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1.5;A_star=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TXT FILES
%fsolve tolerances
options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12);

m=0;
ss=zeros(5000,10);
stable_ss=zeros(2500,4);

k=1;n=4;theta_a=0.5;theta_b=0.5;

l= @(A_star) 1./(1+exp(-(16*A_star-8)));%lambda([ATP])

for b = [0 0.25 0.5 0.75]
  
    fprintf('Running array with a=%.2f, b=%.2f and A*=%.3f at %s.\n',a,b,A_star,datestr(now,'HH:MM:SS'));
   
    syms x1 x2;
    f1 = [l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;...
        l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
    v=[x1,x2];
    jac=jacobian(f1,v); %calculating jacbian with respect to variables x1 & x2

    fhandle=@(X)ode_sigmoid_model(X,a,b,k,A_star,n,theta_a,theta_b);

    for i=0:0.2:2
        for j=0:0.2:2
            m=m+1;
            
            %initial conditions pairing
            X0 = [i,j];
            %using fsolve
            [X,fval,exitflag,output] = fsolve(fhandle,X0,options);

            if (X(1) >=0) && (X(2) >=0)
                %moving row in matrix down by 1 in each loop
                m=m+1;
                
                %rounded calculated steady state values
                x1_ss=round(X(1),2);x2_ss=round(X(2),2);
                
                %subs. in steady state values to jacobian
                sub=subs(jac, [x1 x2], [X(1) X(2)]); %subs. in ss values from original ics
                %calc eigenvlaues
                eigen = eig(sub); %calc eigenvlaues of matrix 'sub'
                %calculate the sign of each eigenvalue
                eigenvalue_1_2=sign(eigen(1));eigenvalue_2_2=sign(eigen(2));
                
                %testing if the steady state is stable or unstable
                if (eigenvalue_1_2 < 0) && (eigenvalue_2_2 < 0)
                    stability = 1;
                else
                    stability = -1;
                end

                %matrix of a, b, energy, ics, ss positions, the stability, lambda and subs2
                ss(m,:) = [a b A_star i j x1_ss x2_ss stability fval(1) fval(2)];
                %ss(~any(ss,2),:) = [];
            else
                disp('negative ss');
            end
        end
    end
end

txtFileName = sprintf('txt-files\\array-b-a=%.0f-energy=%.0f.txt',a*100,A_star*100);
fulltxtFileName=fullfile(txtFileName);
ss(~any(ss,2),:) = [];
fid = fopen(fulltxtFileName,'wt');
for ii = 1:size(ss,1)
    fprintf(fid,'%20.18f\t',ss(ii,:));
    fprintf(fid,'\n');
end

%% FIGURES
fn='Helvetica';wd=16;ht=14;fs_labels=12;fs_axis=12;%fs_title=12;

p=0;

txtFileName = sprintf('txt-files\\array-b-a=%.0f-energy=%.0f.txt',a*100,A_star*100);
ss = importdata(txtFileName);

col_check=ss(:,9);ss_b=ss(col_check < 1e-8,:);
col_check2=ss_b(:,10);ss_c=ss_b(col_check2 < 1e-8,:);
col_check3=ss_c(:,9);ss_d=ss_c(col_check3 > -1e-8,:);
col_check4=ss_d(:,10);ss_1=ss_d(col_check4 > -1e-8,:);

array_fig=figure;
set(gcf,'Units','centimeters','Position',[0 0 wd ht],'PaperUnits','centimeters','PaperSize',[wd ht]);

for b = [0 0.25 0.5 0.75]
    p=p+1;
    
    col_b1=ss_1(:,2);
    ss_2 = ss_1(col_b1 == b,:);
    
    stab_col1=ss_2(:,8);stable_ss1 = ss_2(stab_col1==1,:);unique_stable_ss1 = unique(stable_ss1(:,[1 2 3 6 7]),'rows');
    unstab_col1=ss_2(:,8);unstable_ss1 = ss_2(stab_col1==-1,:);unique_unstable_ss1 = unique(unstable_ss1(:,[1 2 3 6 7]),'rows');
    
% %     disp(unique_stable_ss1);
    
    subplot(2,2,p);
    plot(unique_unstable_ss1(:,4),unique_unstable_ss1(:,5),'ro','MarkerSize',1.8,'MarkerFaceColor','r','MarkerSize',3);hold on;
    plot(unique_stable_ss1(:,4),unique_stable_ss1(:,5),'bo','MarkerSize',1.8,'MarkerFaceColor','b','MarkerSize',3);hold off;
    grid on;box on;ax = gca;ax.XTick = 0:1:3;ax.YTick = 0:1:3;ax.YLim = [0 3];ax.XLim = [0 3];
    ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';
% %     hx=xlabel('$x_1$');hx.Interpreter='latex';hx.FontSize=fs_labels;hx.FontName=fn;%x-axis
% %     hy=ylabel('$x_2$');hy.Interpreter='latex';hy.FontSize=fs_labels;hy.FontName=fn;%y-axis
% %     sub_tit=title(sprintf('b = %.2f',b)); %changing title within for loop for value of ATP used in calculations
% %     sub_tit.FontName=fn;sub_tit.FontWeight='normal';sub_tit.FontSize=12;%fs_axis;

end
    
epsFileName = sprintf('figures\\array-b-a=%.0f-energy=%.0f.eps',a*100,A_star*100);fullFileName =fullfile(epsFileName);print(array_fig,fullFileName,'-depsc');
tiffFileName = sprintf('figures\\array-b-a=%.0f-energy=%.0f.tiff',a*100,A_star*100);fullFileName=fullfile(tiffFileName);print(array_fig,fullFileName,'-dtiff');