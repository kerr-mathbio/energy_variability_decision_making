%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bifurcation diagram time-dependant figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre-setting figure properties
fn='Helvetica';wd=10;ht=9;fs_axis=12;
%fixed parameter values for ODEs
a=1;b=1;k=1;n=4;theta_a=0.5;theta_b=0.5;
%ode45 tolerances
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
%time range to integrate over
tspan=linspace(0,15,10000);

%energy levels to plot diagrams for 
A_star=0.25;%1 stable attractor
% % A_star=0.5;%2 stable attractors
% % A_star=0.8;%3 stable attractors

%lambda(A_star)
l = @(A_star) 1./(1+exp(-(16*A_star-8)));
%ODEs
odes = @(t,x) [l(A_star)*a*x(1)^n./(theta_a^n+x(1)^n)+l(A_star)*b*theta_b^n./(theta_b^n+x(2)^n)-k*x(1);...
                l(A_star)*a*x(2)^n./(theta_a^n+x(2)^n)+l(A_star)*b*theta_b^n./(theta_b^n+x(1)^n)-k*x(2)];
%symbolic variables for protein levels
syms x1 x2;
%ODEs
f_sym = [l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];
% variables for jacobian matrix
v_sym=[x1,x2];
%calculating jacobian with respect to variables x1 & x2
jac=jacobian(f_sym,v_sym);

%creating figure
fig1=figure(1);clf;hold on;
    for i=0:0.3:3%initial conditions on x axis
        for j=0:0.3:3%initial conditions on y axis
            %initial conditions pairing
			ics=[i,j];
			%using ode45
            [t,x_num]=ode45(odes,tspan,ics,options);
			%calculated steady state values
			x1_ss=x_num(2000,1);x2_ss=x_num(2000,2);
			%rounded calculated steady state values
			x1_ss2=round(x1_ss,3);x2_ss2=round(x2_ss,3);

			%subs. in steady state values to jacobian
			sub=subs(jac, [x1 x2], [x1_ss x2_ss]); %subs. in ss values from original ics
			%calc eigenvlaues
			eigen = eig(sub); 
			%calculate the sign of each eigenvalue
			eigenvalue_1=sign(eigen(1));eigenvalue_2=sign(eigen(2));
			
			%testing if the steady state is stable or unstable
			if (eigenvalue_1 < 0) && (eigenvalue_2 < 0)
				stability = 1;
			else 
				stability = -1;
			end

			%x and y axis limits + box and grid for figure
            xlim([0 15]);ylim([0 3]);grid on;box on;
			%plotting time-dependant solution if end steady state is stable
			if stability == 1
				 plot(t,x_num(:,2),'b','LineWidth',0.6);
			else
				point = 0;
			end
            ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';%changing x and y axes properties
            fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];%setting figure size
        end
    end
hold off;
%saving produced figure to output directory with specified name and file extenstion
% % epsFileName = sprintf('figures\\bif_ex_ss3.eps');fullFileName=fullfile(epsFileName);print(fig1,fullFileName,'-depsc');
% % tiffFileName = sprintf('figures\\bif_ex_ss3.tiff');fullFileName2=fullfile(tiffFileName);print(fig1,fullFileName2,'-dtiff');

% % epsFileName = sprintf('figures\\bif_ex_ss2.eps');fullFileName=fullfile(epsFileName);print(fig1,fullFileName,'-depsc');
% % tiffFileName = sprintf('figures\\bif_ex_ss2.tiff');fullFileName2=fullfile(tiffFileName);print(fig1,fullFileName2,'-dtiff');

epsFileName = sprintf('figures\\bif_ex_ss1.eps');fullFileName=fullfile(epsFileName);print(fig1,fullFileName,'-depsc');
tiffFileName = sprintf('figures\\bif_ex_ss1.tiff');fullFileName2=fullfile(tiffFileName);print(fig1,fullFileName2,'-dtiff');
