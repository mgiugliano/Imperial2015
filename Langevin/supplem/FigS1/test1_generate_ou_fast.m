%
% Test suite for generate_ou_fast.m
%
% Sep 2nd 2010 - Michele Giugliano, PhD
%

tau  = 1.;      % eg. units are msec
dt   = tau/80.; % same units of 'tau' 

M    = fix(3 * tau / dt); % transient to be ignored
Npts = M + 1000000.;   % number of points to be generated
R    = 5;           % number of repetition per point

% Desired range of standard deviations
Sd   = [1 2 5 8 10 20 50 80 100 200 500 800 1000 2000 5000 10000]

% Actual estimate of the std and their std of the estimator
Se    = zeros(size(Sd));
sSe   = zeros(size(Sd));
ind   = 1;

for sigma=Sd,
 c1  = 1. - dt / tau;
 c2  = sigma * sqrt(2. * dt / tau);
 tmp = zeros(R,1);

 for h=1:R,
  out = generate_ou_fast(sigma, tau, dt, c1, c2, Npts);
  out = out(M:end);
  tmp(h) = std(out);
 end
  
 Se(ind)  = mean(tmp);
 sSe(ind) = std(tmp);
 ind = ind + 1;
 disp(sprintf('%.0f %% done...', 100. * (sigma-Sd(1))/(Sd(end)-Sd(1))));
end

clf;
figure(1);
hold on;
P = plot(Sd, Sd);
set(P, 'Color', [0 0 0], 'LineWidth', 1);

Q = errorbar(Sd, Se, sSe);
set(Q, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'Color', [0 0 0]);

set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.8 12000], 'YLim', [0.8 12000]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
set(gca, 'XTick', [1 10 100 1000 10000], 'YTick', [1 10 100 1000 10000])
mystr = sprintf('\\sigma_y (estimated over %d points)', Npts - M);
xlabel('\sigma_x', 'FontSize', 20); ylabel(mystr, 'FontSize', 20)
print(gcf, 'panel1.eps', '-loose', '-depsc2');
print(gcf, 'panel1.png', '-loose', '-dpng');



