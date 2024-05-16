clc;
clear;
close all;

%% Burger Equation Parameters
alpha = 1;      % Coefficient of advection term
beta = 10;      % Coefficient of diffusion term

%% Define Domain
L = 2*pi;                      % Length of the spatial domain
x = linspace(-L/2, L/2, 100);  % Spatial discretization

% Initial condition
X0 = cos(x) + sin(x);

% Define the system of ODEs for Burger equation
odefunc = @(t,u) Burger(t,u,x,L,alpha,beta);

% Time discretization
dt = 0.025;                 % Time step
t = 0:dt:9*dt;              % Simulation time

% Solve the system of ODEs Using ODE45
[T,X] = ode45(odefunc, t, X0);

%% Plot the simulation result

% Create a meshgrid for 3D plotting
[Tm,Xm] = meshgrid(t,x);
% Plot the surface
surf(Xm, Tm, X.');
xlabel('x');
ylabel('t');
zlabel('u');
title('Burger Equation Simulation');
