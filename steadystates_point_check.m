%%%code in this script is produced with comments explaining what the line of code next to or below it does

%parameter values for ODEs
k=1;
theta_a=0.5;theta_b=0.5;
n=4;
a=1;b=1.5;

%lambda, if linear function is chosen
% % l= @(A_star) A_star;
%lambda, if sigmoid function is chosen
l= @(A_star) 1./(1+exp(-(16*A_star-8)));
%energy value
A_star=0.9;

%calculated steady state
steadystate=[0.8593764644e-1, 2.909085952];

%symbolic variables for protein levels
syms x1 x2;
%ODEs
ODEs=[l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2];

%sub steady state back into ODEs for check
sub_into_odes = subs(ODEs,[x1 x2],steadystate);
%converting to floating point number
sub_into_odes2=vpa(sub_into_odes);

%variables for jacobian
variables=[x1,x2];
%calculating jacbian with respect to variables x1 & x2
jac=jacobian(ODEs,variables); 
%substitute steady state into Jacobian matrix
sub=subs(jac, [x1 x2], steadystate);
%calculate eigenvalues
eigenvalues = eig(sub);
%sign of both eigenvalues
eigen_sign=([sign(eigenvalues(1)),sign(eigenvalues(2))]);



%displaying sub_into_odes
disp('Steady state evaluated in ODEs =');
disp(sub_into_odes2);

%display sign of both eigenvalues in command window
disp('Eigenvalues sign =');
disp(eigen_sign);
