function out = generate_ou_fast(sigma, tau, dt, c1, c2, Npts)
%
% y = generate_ou_fast(sigma_x, tau_x, dt, speed_up, Npts)
%
% Generates a realisation of a discrete-time Ornstein-Uhlenbeck
% process, employing NOT the Gillespie, DT (1996) exact algorithm
% but its Taylor-expanded Euler approximation.
%
% At the steady-state (reached after a time of the order ot 'tau'),
% the process will have zero mean, variance 'sigma^2' and covariance
% exponentially decaying with a single time-constant 'tau'
%
% y: [Npts x 1] output vector 
%
% sigma:   steady-state desired standard deviation of y
% tau:     steady-state autocorrelation time-length of y
% dt:      iteration time-step [same units of 'tau']
% c1:      pre-calculated coefficient '1 - dt / tau'
% c2:      pre-calculated coefficient 'sigma * sqrt(2 * dt / tau)'
% Npts:    number of points to be generated
%
%
% Sep 2nd 2010 - Michele Giugliano, PhD
%

out = zeros(Npts, 1);     % Memory pre-allocation (for speed purposes)

y = 0.;                   % Initial condition

for k=1:Npts,
 out(k) = y; 
 y = c1 * y + c2 * randn;
end

