%declaring a function that accepts inputs X, a, b, k, energy, n, theta_a, theta_b and returns outputs in F
function F = ode_sigmoid_model(X,a,b,k,A_star,n,theta_a,theta_b)

%setting variable names for protein levels in ODEs
x1=X(1);
x2=X(2);

%lambda coefficients & function
l = @(A_star) 1./(1+exp(-(16*A_star-8)));

%ODEs for x1 and x2 with time
x1dot=l(A_star)*a*x1^n./(theta_a^n+x1^n)+l(A_star)*b*theta_b^n./(theta_b^n+x2^n)-k*x1;
x2dot=l(A_star)*a*x2^n./(theta_a^n+x2^n)+l(A_star)*b*theta_b^n./(theta_b^n+x1^n)-k*x2;

%ODE system of equations
F=[x1dot;x2dot];