% Quick simulation to test cr3bp.m and jacobi.m
clear,clf,clc
%% Problem Setup
mu = 1/82.3; % Earth-Moon
T  = 2*pi;   % nondim. period of primaries

[L1,gamma,xL] = findGamma(1,mu);

XL = [xL(1); 0; 0; 0];

%% Eigendecomposition of Planar Phase-Space

% L1 eigenbasis
mubar = mu/abs(xL-1+mu)^3 + (1-mu)/abs(xL+mu)^3;

a = 1 + 2*mubar;
b = mubar - 1;

lam = sqrt( 0.5*(mubar-2+sqrt(9*mubar^2-8*mubar)));
nu  = sqrt(-0.5*(mubar-2-sqrt(9*mubar^2-8*mubar)));

sigma = 2*lam/(lam^2 + b);
tau   =-(nu^2 + a)/(2*nu); 

% basis vectors
u1 = [1; -sigma;  lam; -lam*sigma];
u2 = [1;  sigma; -lam; -lam*sigma];
u  = [1;    0  ;   0 ;   nu*tau  ];
v  = [0;   tau ; -nu ;     0     ];
L1eig = [u1 u2 u v];
%% Initial Conditions 1
displacement = 1e-5;
X0 = XL + displacement*u;

% planar Lyapunov orbit period
tfinal = 2*pi/nu; 

tspan = 0:0.001*T:tfinal;

%% Numerical Integration
tic
opts = odeset('RelTol',1e-13,'AbsTol',1e-16);
[t,X] = ode113(@(t,x) cr3bp(t,x,mu),tspan,X0,opts);
toc

x  = X(:,1);
y  = X(:,2);
vx = X(:,3);
vy = X(:,4);

%% Plotting
cr3bp_plot(t,X,mu,Lpt,1)

disp(['Numerical Error (Energy Conservation): ' num2str(std(E))])