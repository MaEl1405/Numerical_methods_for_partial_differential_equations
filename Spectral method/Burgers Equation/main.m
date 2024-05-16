clc;clear; close all;

%% Burger Equation parameters
alpha = 1  ;             % Diffusion Constant
beta  = 10 ;

%% Define Spatial Domain
L  = 2*pi;                      % Length of Domain
N  = 100;                       % Number of Discritization POint
dx = L / N;     
x  = -L / 2:dx:L / 2 - dx;      % X domain


%% Define Descrite Wavenumber
kappa = (2 * pi / L) * (-N / 2 : N / 2 - 1);
kappa = fftshift(kappa');

%% Initial Condition 
u0 = cos(x)+sin(x);

%% Solve Converted PDE to n-dimensions ODE
dt = 0.025;
t = 0:dt:9*dt;
[t, u] = ode45(@(t, u)rhsBurgers(t, u, kappa, alpha,beta,A,w,x), t, u0); % Integrate up to t=10 exactly

%% Plot Results
figure
surf(x, t, real(u))
xlabel('x',FontSize=20)
ylabel('t(s)',FontSize=20)
zlabel('u(x, t)',FontSize=20)
title('Solution of the 1D Burgers'' Equation')

grid on
