%
% POTASSIUM CURRENTS
%

syms alpha_n beta_n N;
syms mean_true variance_true;
syms mean_Fox variance_Fox;


%---- FOX ---
syms mn vn mn4 vn4;

mn = alpha_n / (alpha_n + beta_n);
vn = alpha_n*beta_n / (N * (alpha_n + beta_n)^2);

mn4 = mn^4 + 6 * mn^2 * vn + 3 * vn^2;
vn4 = vn * (16 * mn^6 + 168*mn^4*vn + 384 * mn^2 * vn^2 + 97 * vn^3);

mean_Fox = mn4;
variance_Fox = vn4;
%------


%--- TRUTH
syms nbar v1 v2 v3 v4;
nbar = (alpha_n / (alpha_n + beta_n));
v1 = 4 * nbar^7 * (1 - nbar) / N;
v2 = 6 * nbar^6 * (1 - nbar)^2 / N;
v3 = 4 * nbar^5 * (1 - nbar)^3 / N;
v4 = nbar^4 * (1 - nbar)^4 / N;
mean_true = (alpha_n / (alpha_n + beta_n))^4;
variance_true = v1 + v2 + v3 + v4;

%----


% Desired range of Voltages
Vd   =  -60:2.5:50;

m_true = zeros(size(Vd));
s_true = zeros(size(Vd));
m_Fox  = zeros(size(Vd));
s_Fox  = zeros(size(Vd));

N = 100;

ind = 1;
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
print(gcf, 'panel1.eps', '-loose', '-depsc2');
print(gcf, 'panel1.png', '-loose', '-dpng');



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
print(gcf, 'panel2.eps', '-loose', '-depsc2');
print(gcf, 'panel2.png', '-loose', '-dpng');



