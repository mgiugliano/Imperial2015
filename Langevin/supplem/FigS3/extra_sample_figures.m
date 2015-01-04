load sample_data.mat; % out @ V = 50 mV
dt= 0.001;

xout = out.^4;


 [c lags] = xcov(xout, fix(15. / dt), 'coeff');  N = length(c);
  X = dt*lags(1+(N-1)/2:end)';
  Y = c(1+(N-1)/2:end);
 
  s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0],...
               'Upper',[Inf],...
               'Startpoint',[1]);
f = fittype('exp(-x/a)','options',s);

  [model,gof2] = fit(X,Y,f);
figure(25);
cla;

hold on;
Q = plot(X, exp(-X/model.a));
%Q = plot(model,'m'); 
set(Q, 'Color', [0.8 0 0], 'LineStyle', '-', 'LineWidth', 2);
P = plot(X(1:100:end), Y(1:100:end), 'o');
set(P, 'MarkerFaceCOlor', 'none', 'MarkerEdgeCOlor', [0 0 0])
hold off;

set(gca, 'XLim', [0 15]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
xlabel('\Delta (msec)', 'FontSize', 20); ylabel('normalized covariance', 'FontSize', 20)
print(gcf, 'panel_extra2.eps', '-loose', '-depsc2');
print(gcf, 'panel_extra2.png', '-loose', '-dpng');


figure(26); cla;
[NNN XXX] = hist(xout,1000);
dX = XXX(2)-XXX(1);
NNN = NNN / (length(out) * dX);
hold on;
%B = bar(XXX, NNN, 1); set(B, 'EdgeColor', [0 0 0], 'FaceCOlor', [0 0 0]);
B = plot(XXX(1:10:end), NNN(1:10:end), 'o'); set(B, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]) 
Q = plot(XXX, (1./sqrt(2*pi*std(xout)^2)*exp(-(XXX - mean(xout)).^2./(2*std(xout)^2)))); 
set(Q, 'Color', [0.8 0 0], 'LineStyle', '-', 'LineWidth', 2);
hold off;

%set(gca, 'XLim', [0.9 1.04]);
set(gca, 'FontName', 'Arial', 'FontSize', 15, 'XGrid', 'on', 'YGrid', 'on', 'box', 'on');
xlabel('u', 'FontSize', 20); %ylabel('normalized covariance', 'FontSize', 20)
print(gcf, 'panel_extra.eps', '-loose', '-depsc2');
print(gcf, 'panel_extra.png', '-loose', '-dpng');
