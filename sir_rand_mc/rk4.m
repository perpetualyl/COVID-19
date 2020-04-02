function[Y]=rk4(y, yprime, t0, dt, N);
% Y = rk4(y, yprime, t0, dt, N, history=true)
%
% Solves the ODE system by taking N timesteps of size dt using the fourth-order
% Runge-Kutta method. 
%
% If the optional input history is false, then only the terminal time state is
% output; otherwise the entire time history is saved.

if nargin < 6
  history = true;
end

t = t0;

if history
  Y = zeros([numel(y) N+1]);
  Y(:,1) = y;

  for k = 1:N
    k1 = yprime(t, Y(:,k));
    k2 = yprime(t+dt/2, Y(:,k) + dt/2*k1);
    k3 = yprime(t+dt/2, Y(:,k) + dt/2*k2);
    k4 = yprime(t+dt, Y(:,k) + dt*k3);

    Y(:,k+1)=Y(:,k) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    t = t + dt;
  end

else
  Y = y;
  for k = 1:N
    k1 = yprime(t, Y);
    k2 = yprime(t+dt/2, Y + dt/2*k1);
    k3 = yprime(t+dt/2, Y + dt/2*k2);
    k4 = yprime(t+dt, Y + dt*k3);

    Y = Y + dt/6*(k1 + 2*k2 + 2*k3 + k4);
  end
end
