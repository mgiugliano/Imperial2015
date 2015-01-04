%
% Test suite for generate_Fox_fast.m
%
% Sep 2nd 2010 - Michele Giugliano, PhD
%


s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...
               'Upper',[Inf],...
               'Startpoint',[1]);
f = fittype('exp(-x/a)','options',s);


dt   = 0.001; % same units of 1/alpha or 1/beta 

M    = fix(3 * 30. / dt); % transient to be ignored
Npts = M + 10000000.;   % number of points to be generated
R    = 5;           % number of repetition per point

% Desired range of Voltages
Vd   =  -60:5:50;

% Actual estimate of the std and their std of the estimator
tUe    = zeros(size(Vd));
stUe   = zeros(size(Vd));

theory_tUe = zeros(size(Vd));

ind   = 1;

for V=Vd,
% CHOSE ONE SET OF ALPHA, BETA FROM BELOW ---------------------------------
%
%
%if (V==-40), alpha_m = 1; else alpha_m = -0.1 * (V+40.)/(exp(-0.1*(V+40.)) -1.); end;
%beta_m  =   4. * exp(-(V+65.)/18.);
%alpha_h =   0.07 * exp(-(V+65.)/20.);
%beta_h  =   1. / (exp(-0.1 * (V+35.)) + 1.);
%if (V==-55), alpha_n = 0.1; else alpha_n = -0.01 * (V+55.)/ ( exp(-0.1*(V+55.)) - 1. ); end
%beta_n  =   0.125 * exp(-(V+65.)/80.);

if (V==-55), alpha_n = 0.1; else alpha_n = -0.01 * (V+55.)/ ( exp(-0.1*(V+55.)) - 1. ); end
beta_n  =   0.125 * exp(-(V+65.)/80.);
alpha = alpha_n;
beta  = beta_n;
%--------------------------------------------------------------------------
 
    
 tmp1 = zeros(R,1);
 Nchan = 100.;
 
 for h=1:R,
  out = generate_Fox_fast(V, dt, Nchan, alpha, beta, Npts);
  out = out(M:end);

  [c lags] = xcov(out, fix(15. / dt), 'coeff');  N = length(c);
  X = dt*lags(1+(N-1)/2:end)';
  Y = c(1+(N-1)/2:end);
 
  [model,gof2] = fit(X,Y,f);
%   figure(1); clf;
%   plot(X, Y)
%   hold on; plot(model,'m'); hold off;
%   pause;
  tmp1(h) = model.a;
 end
  
 tUe(ind)   = mean(tmp1);
 stUe(ind)  = std(tmp1);
 theory_tUe(ind) = 1./(alpha+beta);
 
 ind = ind + 1;
 disp(sprintf('%.0f %% done...', 100. * (V-Vd(1))/(Vd(end)-Vd(1))));
end

clf;
figure(1);
hold on;
P = plot(Vd, theory_tUe);
set(P, 'Color', [0 0 0], 'LineWidth', 1);

Q = errorbar(Vd, tUe, stUe);
set(Q, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'Color', [0 0 0]);

set(gca, 'XLim', [-65 55], 'XTick', [-60:10:50]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
xlabel('V', 'FontSize', 20); ylabel('\tau_u', 'FontSize', 20)
print(gcf, 'panel3.eps', '-loose', '-depsc2');
print(gcf, 'panel3.png', '-loose', '-dpng');


figure(25);
cla;
hold on;
Q = plot(model,'m'); 
set(Q, 'Color', [0.6 0.6 0.6], 'LineStyle', '-', 'LineWidth', 2);
P = plot(X(1:100:end), Y(1:100:end), 'o');
set(P, 'MarkerFaceCOlor', 'none', 'MarkerEdgeCOlor', [0 0 0])
hold off;

set(gca, 'XLim', [0 15]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
xlabel('\Delta (msec)', 'FontSize', 20); ylabel('normalized covariance', 'FontSize', 20)
print(gcf, 'panel4.eps', '-loose', '-depsc2');
print(gcf, 'panel4.png', '-loose', '-dpng');


figure(26); cla;
[NNN XXX] = hist(out,1000);
dX = XXX(2)-XXX(1);
NNN = NNN / (length(out) * dX);
hold on;
%B = bar(XXX, NNN, 1); set(B, 'EdgeColor', [0 0 0], 'FaceCOlor', [0 0 0]);
B = plot(XXX(1:10:end), NNN(1:10:end), 'o'); set(B, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]) 
Q = plot(XXX, (1./sqrt(2*pi*std(out)^2)*exp(-(XXX - mean(out)).^2./(2*std(out)^2)))); 
set(Q, 'Color', [0.8 0 0], 'LineStyle', '-', 'LineWidth', 2);
hold off;

set(gca, 'XLim', [0.9 1.04]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
xlabel('u', 'FontSize', 20); %ylabel('normalized covariance', 'FontSize', 20)
print(gcf, 'panel5.eps', '-loose', '-depsc2');
print(gcf, 'panel5.png', '-loose', '-dpng');
