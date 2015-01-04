%
% SODIUM CURRENTS
%

syms alpha_m beta_m alpha_h beta_h N;
syms mean_true variance_true;
syms mean_Fox variance_Fox;


%---- FOX ---
syms mh vh mm vm mm3 vm3;

mm = alpha_m / (alpha_m + beta_m);
vm = alpha_m*beta_m / (N * (alpha_m + beta_m)^2);

mh = alpha_h / (alpha_h + beta_h);
vh = alpha_h*beta_h / (N * (alpha_h + beta_h)^2);

%mm3= alpha_m^2 * (alpha_m + 3 * beta_m / N) / (alpha_h + beta_h)^3;
%vm3= 3*alpha_m^3*beta_m*(3*alpha_m^2 + 12 * alpha_m*beta_m/N + 5*beta_m^2/(N^2))/(N*(alpha_m+beta_m)^6);
mm3 = mm * (mm^2 + 3 * vm);
vm3 = 3 * vm * (3 * mm^4 + 12 * mm^2 * vm + 5 * vm^2);
mean_Fox = mm3 * mh;
variance_Fox = vm3 * vh + vm3 * mh^2 + vh * mm3^2;
%------


%--- TRUTH
syms mbar hbar v1 v2 v3 v4 v5 v6 v7;
mbar = (alpha_m / (alpha_m + beta_m));
hbar = (alpha_h / (alpha_h + beta_h));
v1   = mbar^6 * hbar * (1 - hbar) / N;
v2   = 3 * mbar^5 * hbar^2 * (1 - mbar) / N;
v3   = 3 * mbar^4 * hbar^2 * (1 - mbar)^2 / N;
v4   = mbar^3 * hbar^2 * (1 - mbar)^3 / N;
v5   = 3 * mbar^5 * hbar * (1 - mbar) * (1 - hbar) / N;
v6   = 3 * mbar^4 * hbar * (1 - mbar)^2 * (1 - hbar) / N;
v7   =  mbar^3 * hbar * (1 - mbar)^3 * (1 - hbar) / N;
mean_true = (alpha_m / (alpha_m + beta_m))^3 * (alpha_h / (alpha_h + beta_h));
variance_true = v1 + v2 + v3 + v4 + v5 + v6 + v7;

%----




% Desired range of Voltages
Vd   =  -60:2.5:50;

m_true = zeros(size(Vd));
s_true = zeros(size(Vd));
m_Fox  = zeros(size(Vd));
s_Fox  = zeros(size(Vd));

N = 1200;

ind = 1;
for V=Vd,

%--------------------------------------------------------------------------
if (V==-40), alpha_m = 1; else alpha_m = -0.1 * (V+40.)/(exp(-0.1*(V+40.)) -1.); end;
beta_m  =   4. * exp(-(V+65.)/18.);
alpha_h =   0.07 * exp(-(V+65.)/20.);
beta_h  =   1. / (exp(-0.1 * (V+35.)) + 1.);
%--------------------------------------------------------------------------
m_true(ind) = eval(mean_true);
s_true(ind) = sqrt(eval(variance_true));
m_Fox(ind)  = eval(mean_Fox);
s_Fox(ind)  = sqrt(eval(variance_Fox));

ind = ind + 1;
end

figure(1); %clf;
hold on;
P = plot(Vd, m_true);
set(P, 'Color', [0 0 0], 'LineWidth', 1);
Q = plot(Vd, m_Fox, '--');
set(Q, 'Color', [1 0 0], 'LineWidth', 1);
hold off;

set(gca, 'XLim', [-65 55], 'XTick', [-60:10:50]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
xlabel('V  (mV)', 'FontSize', 20); ylabel('mean', 'FontSize', 20)
L = legend([P Q], 'Theory', 'Fox (1997)'); legend boxoff; set(L, 'FontSize', 20)
print(gcf, 'panel3.eps', '-loose', '-depsc2');
print(gcf, 'panel3.png', '-loose', '-dpng');



figure(2); %clf;
hold on;
P = plot(Vd, s_true);
set(P, 'Color', [0 0 0], 'LineWidth', 1);
Q = plot(Vd, s_Fox, '--');
set(Q, 'Color', [1 0 0], 'LineWidth', 1);
hold off;

set(gca, 'XLim', [-65 55], 'XTick', [-60:10:50]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
xlabel('V  (mV)', 'FontSize', 20); ylabel('std', 'FontSize', 20)
L = legend([P Q], 'Theory', 'Fox (1997)'); legend boxoff; set(L, 'FontSize', 20)
print(gcf, 'panel4.eps', '-loose', '-depsc2');
print(gcf, 'panel4.png', '-loose', '-dpng');



