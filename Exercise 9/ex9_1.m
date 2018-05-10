%% Initialize varables
clearvars
close all

% Node link incidence matrix
B = [1 1 0 0;
     -1 0 1 0;
     0 -1 0 1;
     0 0 -1 -1];

% Exogenous inflows
lambda = zeros(4,1);
lambda(1) = 1;

% Exogenous outflows
mu = zeros(4,1);
mu(4) = 1;

% Net exogenous flow
nu = lambda - mu;

%% a) Social optimum flow calculation
cvx_begin
    variable f(4)
    minimize f(1)*2*f(1) + f(2)*2 + f(3)*3 + f(4)*3*f(4)
    subject to
        B*f == nu
        f >= 0
cvx_end
 
fprintf('Social optimum flow: \n%8s %6s %6s %6s %6s \n','Link','1','2','3','4')
fprintf('%8s %6.2f %6.2f %6.2f %6.2f \n','Flow',f')

%% b) User optimum flow calculation (Wardrop equilibrium)
cvx_begin
    variable f(4)
    minimize pow_p(f(1),2) + 2*f(2) + 3*f(3) + 3/2*pow_p(f(4),2)
    subject to
        B*f == nu
        f >= 0
cvx_end

fprintf('User optimum flow: \n%8s %6s %6s %6s %6s \n','Link','1','2','3','4')
fprintf('%8s %6.2f %6.2f %6.2f %6.2f \n','Flow',f')

%% c) Tolls
w = zeros(4,1);
w(1) = 1;
w(4) = 3/2;

% User optimum
cvx_begin
    variable f(4)
    minimize pow_p(f(1),2) + w(1)*f(1) ...
        + 2*pow_p(f(2),1) + w(2)*f(2) ...
        + 3*pow_p(f(3),1) + w(3)*f(3) ...
        + 3/2*pow_p(f(4),2) + w(4)*f(4)
    subject to
        B*f == nu
        f >= 0
cvx_end

fprintf('User optimum flow: \n%8s %6s %6s %6s %6s \n','Link','1','2','3','4')
fprintf('%8s %6.2f %6.2f %6.2f %6.2f \n','Flow',f')




