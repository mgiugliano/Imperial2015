function out = generate_Fox_fast(V, dt, N, alpha, beta, Npts)
%
% y = generate_Fox_fast(V, dt, N, alpha, beta, Npts)
%
% Generates a realisation of a discrete-time Fox's 
% process, employing NOT the Gillespie, DT (1996) exact algorithm
% but its Taylor-expanded Euler approximation.
%
% At the steady-state (reached after a time of the order ot '1/(alpha+beta)'),
% the process will have mean (alpha/(alpha+beta)), variance 'sigma^2' (see 
% the text) and covariance exponentially decaying with a single time-constant
% '1/(alpha+beta)'.
%
% y: [Npts x 1] output vector 
%
% V:       clamped membrane potential
% dt:      iteration time-step [same units of 'tau']
% Npts:    number of points to be generated
%
%
% Sep 2nd 2010 - Michele Giugliano, PhD
%

c1 = 1 - (alpha + beta) * dt;
c2 = alpha * dt;
c3 = sqrt(2 * alpha * beta * dt / (N * (alpha + beta)));

out = zeros(Npts, 1);     % Memory pre-allocation (for speed purposes)

y = 0.;                   % Initial condition

for k=1:Npts,
 out(k) = y; 
 y = c1 * y + c2 + c3 * randn;
end

