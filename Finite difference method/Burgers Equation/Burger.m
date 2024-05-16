function dudt = Burger(t, u, x, L, alpha, beta)
% BURGER computes the time derivative of the state vector u for the Burger equation.
%   dudt = BURGER(t, u, x, L, alpha, beta) computes the time derivative of the
%   state vector u for the Burger equation at time t and spatial variable x.
%
%   Inputs:
%   t      : Current time
%   u      : State vector representing the dependent variable (e.g., velocity) at spatial locations x
%   x      : Spatial variable
%   L      : Length of the spatial domain
%   alpha  : Coefficient of the advection term
%   beta   : Coefficient of the diffusion term
%
%   Output:
%   dudt   : Time derivative of u according to the Burger equation

    % Length of the state vector u
    N = length(u);
    
    % Identity matrix of size N
    I = eye(N);
    
    % Matrix representing the first derivative with periodic boundary conditions
    A1 = circshift(I, -1) - circshift(I, 1);
    
    % Matrix representing the second derivative with periodic boundary conditions
    A2 = circshift(I, 1) + circshift(I, -1) - 2 * I;

    % Spatial grid spacing
    dx = L / N;
    
    % First spatial derivative of u
    ux = (A1 * u) / (2 * dx);
    
    % Second spatial derivative of u
    uxx = (A2 / (dx^2)) * u;
    
    % Time derivative of u computed using the Burger equation
    dudt = alpha * u .* ux + beta * uxx;
end
