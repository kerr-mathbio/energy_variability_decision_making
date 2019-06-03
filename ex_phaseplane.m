%%%code in this script is produced with comments explaining what the line of code next to or below it does

%%%code for 1,2,3 or 4 stable steady states needs to be selected and the remaining commented out

%pre-setting figure properties
fn='Helvetica';wd=10;ht=9;
%setting time range to integrate over
tspan=linspace(0,100,50000);
%parameter values for ODEs that are being fixed
k=1;n=4;theta_a=0.5;theta_b=0.5;
%ode45 tolerances
ode_options=odeset('RelTol',1e-8,'AbsTol',1e-8);

%creating figure
fig1=figure;
%setting matrix row value to zero
matrix_row=0;

%setting parameter sets to get 1-4 stable steady states (sss)
a=1;b=1;A_star=0.2;%1sss
% a=1;b=1.75;A_star=0.5;%2sss
% a=1;b=1;A_star=0.8;%3sss
% a=1;b=0;A_rrstar=1;%4sss

l= @(A_star) 1./(1+exp(-(16*A_star-8)));

%pre-setting matrix size to speed up computations
M1=zeros(150,5);
%ODEs
f = @(t,x) [l(A_star)*a*x(1)^n./(theta_a^n+x(1)^n)+l(A_star)*b*theta_b^n./(theta_b^n+x(2)^n)-k*x(1);...
    l(A_star)*a*x(2)^n./(theta_a^n+x(2)^n)+l(A_star)*b*theta_b^n./(theta_b^n+x(1)^n)-k*x(2)];
%symbolic variables for protein levels
syms x1 x2;
%ODEs
f_sym = [l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
%variables for jacobian matrix
v_sym=[x1,x2];
%calculating jacbian with respect to variables x1 & x2
jac=jacobian(f_sym,v_sym);

%initial conditions on x axis
for i=0:0.3:3
    %display where the computation is up to in command window - good for long computations to see where you are up to
    fprintf('Running subplot with i=%.2f.\n',i);
    %initial conditions on x axis
    for j=0:0.3:3
        %increasing the value by 1 in each new loop - moves the row
        %values up by one each loop
        matrix_row=matrix_row+1;
        %using ode45
        [t,x_num]=ode45(f,tspan,[i,j],ode_options); %solving ODEs with ics
        %calculated steady state values
        x1_ss=x_num(5000,1);x2_ss=x_num(5000,2);
        %rounded calculated steady state values
        x1_ss_b=round(x1_ss,3);x2_ss_b=round(x2_ss,3); %rounding steady state(ss) position to see unique ss
        
        %subs. in steady state values to jacobian
        sub=subs(jac, [x1 x2], [x1_ss x2_ss]); %subs. in ss values from original ics
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
        
        %matrix of ics, ss positions and the stability
        M1(matrix_row,:) = [i j x1_ss_b x2_ss_b stability];
        
        %axes limits & box around figure
        grid on;hold on;box on;
        
        %plotting solution trajectory if stability is stable
        if stability == 1
                        plot(x_num(:,1),x_num(:,2),'Color',[0.1,0.9,1]);%1ss
% %                         plot(x_num(:,1),x_num(:,2),'Color',[0,0.5,1]);%2ss
% %                         plot(x_num(:,1),x_num(:,2),'Color',[0.4,0.1,1]);%3ss
% %             plot(x_num(:,1),x_num(:,2),'m');%4ss
        else
            point = 0;
        end
    end
end

%selecting a column in M1 matrix
col_stable=M1(:,5);
%new sub-matrix M1_b is a submatrix of M1 with stability value = 1 (stable stedy states)
M1_b=M1(col_stable==1,:);
%extracting unique stable steady states
M2 = unique(M1_b(:,[3 4]),'rows');

%plotting stable steady states
plot(M2(:,1),M2(:,2),'o','MarkerEdgeColor','black', 'MarkerFaceColor',[0.49 1 0.63],...
    'MarkerSize',6);
%removing axis labels
set(gca,'YTickLabel',[],'XTickLabel',[]);
%setting figure size
set(gcf,'Units','centimeters','Position',[0 0 wd ht],'PaperUnits','centimeters','PaperSize',[wd ht]);
%setting axis limits
xlim([0 3]);ylim([0 3]);hold off;

%saving produced figure to output directory with specified name and file extenstion

epsFileName1 = 'figures\\1_steady_state.eps';fullFileName =fullfile(epsFileName1);print(fig1,fullFileName,'-depsc');%1sss
tiffFileName1 = 'figures\\1_steady_state.tiff';fullFileName =fullfile(tiffFileName1);print(fig1,fullFileName,'-dtiff');%1sss

% % epsFileName1 = 'figures\\2_steady_states.eps';fullFileName =fullfile(epsFileName1);print(fig1,fullFileName,'-depsc');%2sss
% % tiffFileName1 = 'figures\\2_steady_state.tiff';fullFileName =fullfile(tiffFileName1);print(fig1,fullFileName,'-dtiff');%2sss

% % epsFileName1 = 'figures\\3_steady_states.eps';fullFileName =fullfile(epsFileName1);print(fig1,fullFileName,'-depsc');%3sss
% % tiffFileName1 = 'figures\\3_steady_state.tiff';fullFileName =fullfile(tiffFileName1);print(fig1,fullFileName,'-dtiff');%3sss

% % epsFileName1 = 'figures\\4_steady_states.eps';fullFileName =fullfile(epsFileName1);print(fig1,fullFileName,'-depsc');%4sss
% % tiffFileName1 = 'figures\\4_steady_state.tiff';fullFileName =fullfile(tiffFileName1);print(fig1,fullFileName,'-dtiff');%4sss