clear all, close all, clc

addpath("./utils")
addpath(genpath("./spotless/"))

% Dimension of the problem
d = 3;

% Claim the variable
x = msspoly('x',d);
% Inequality constraints
g = [20 - sum(x.^2)]; % Assume bounded

g = [g;1 - sum(x.^4)]; % TV Screen Example


problem.vars = x;
problem.inequality = g;
% Relaxation Order
kappa = 2;

%% Solve GRCC problem
P = eye(d);  % The projection matrix, identity if the problem is full dimension 
Q = eye(d);  % The shape matrix of the ellipsoid, identity if we want enclosing ball
sol = GRCC(problem,kappa,P,Q,d);

grcc_rad = sol.upper_bound;
real_rad = sqrt(3*(nthroot(1/3,4))^2);

fprintf("\n GRCC rad: %.2f\n Real rad: %.2f \n", grcc_rad, real_rad);