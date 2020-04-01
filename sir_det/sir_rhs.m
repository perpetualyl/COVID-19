function[rhs] = sir_rhs(u, alpha)
% rhs = sir_rhs(u, alpha)
% 
% Evaluates right-hand side of the normalized SIR model. This is the ODE system
%
%   S' = -alpha*S*I
%   I' = alpha*S*I - I
%   R' = I
%
% The input u is a 3 x 1 vector with the components u = [S; I; R].
%
% The un-normalized SIR model has two parameters gamma and beta. This
% normalized model has parameter alpha satisyfing, alpha = beta/gamma, and has
% time rescaled as time(normalized) = gamma*time(un-normalized).

rhs = [-alpha*u(1,:).*u(2,:); ...
       alpha*u(1,:).*u(2,:) - u(2,:); ...
       u(2,:)];
