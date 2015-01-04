%
% Test suite for generate_ou_fast.m
%
% Sep 2nd 2010 - Michele Giugliano, PhD
%


sigma = 50.;
R    = 5;           % number of repetition per point

s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...
               'Upper',[Inf],...
               'Startpoint',[1]);
f = fittype('exp(-x/a)','options',s);

% Desired range of autocorrelation times 
Td    = [0.01 0.02 0.05 0.08 0.1 0.2 0.5 0.8 1 2 5 8 10 20 50 80]

% Actual estimate of the autocorrelation time and the std of the estimator
Te    = zeros(size(Sd));
sTe   = zeros(size(Sd));
ind   = 1;


for tau=Td,
 dt   = tau/80.; % same units of 'tau' 
 M    = fix(3 * tau / dt); % transient to be ignored
 Npts = M + 1000000.;   % number of points to be generated
    
 c1  = 1. - dt / tau;
 c2  = sigma * sqrt(2. * dt / tau);
 tmp = zeros(R,1);

 for h=1:R,
  out = generate_ou_fast(sigma, tau, dt, c1, c2, Npts);
  out = out(M:end);
  
  [c lags] = xcov(out, 3*M, 'coeff');
  N = length(c);
  X = dt*lags(1+(N-1)/2:end)';
  Y = c(1+(N-1)/2:end);
 
  [model,gof2] = fit(X,Y,f);
  %figure(1); clf;
  %plot(X, Y)
  %hold on; plot(model,'m'); hold off;
  %pause;
  tmp(h) = model.a;
 end
  
 Te(ind)  = mean(tmp);
 sTe(ind) = std(tmp);
 ind = ind + 1;
 disp(sprintf('%.0f %% done...', 100. * (tau-Td(1))/(Td(end)-Td(1))));
end

clf;
figure(1);
hold on;
P = plot(Td, Td);
set(P, 'Color', [0 0 0], 'LineWidth', 1);

Q = errorbar(Td, Te, sTe);
set(Q, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'Color', [0 0 0]);

set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [0.008 120], 'YLim', [0.008 120]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
set(gca, 'XTick', [0.01 0.1 1 10 100], 'YTick', [0.01 0.1 1 10 100])

mystr = sprintf('\\tau_y (estimated over %d points)', Npts - M);
xlabel('\tau_x', 'FontSize', 20); ylabel(mystr, 'FontSize', 20)
print(gcf, 'panel2.eps', '-loose', '-depsc2');
print(gcf, 'panel2.png', '-loose', '-dpng');



