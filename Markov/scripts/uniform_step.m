%
% uniform_step.m
%
%
% (c) 2015 - Michele Giugliano, michele.giugliano@uantwerpen.be
% http://www.uantwerpen.be/michele-giugliano
%
% Reference: Linaro D, Giugliano M (2014) Markov models of ion channels. 
% Springer Encyclopedia of Computational Neuroscience
% doi:10.1007/978-1-4614-7320-6_131-1
%

%clear all;      % Clear all varibles from the memory
%close all;      % Close all files and figures
%clc;            % Clear the command window

% Comment the following line, if you want to launch 'benchmark_me.m'
N = 100;       % Number of independent and identical channels

Nstep = 5000;    % Simulation time steps (exclusively used to plot traces)
Nadd  = 10000;   % Additional simulation time steps, until steady-state reached!

dt= 0.0001;      % ms - simulation time step 

a = 1;          % \alpha - [ms^-1] - values of the parameters
b = 7;          % \beta  - [ms^-1] - values of the parameters

STATES = {'A', 'B', 'C', 'D', 'E', 'F'};    % Aesthetics: name of states

% Equivalent initial condition 
xi = 3 * ones(N,1); % equivalent to P = [0 0 1 0 0]';, xi \in {1,2,3,4,5} 

% Generator matrix of A <-> B <-> C <-> D <-> E, like "Kd" ion channels
Q = [-4*a     b         0      0     0; 
      4*a -(3*a+b)     2*b     0     0;
       0     3*a   -(2*a+2*b) 3*b    0;
       0      0        2*a -(a+3*b) 4*b;
       0      0         0      a   -4*b];
      
temp = eye(5) + Q*dt;        % Useful temporary matrix

out = zeros(Nstep, 6);      % Output data structure, for plotting 

for k=1:Nstep,
 for i=1:N,
     xi(i) = RV(temp(:,xi(i))');
 end
 out(k,1)   = k*dt;
 [NN XX] = histc(xi, [1 2 3 4 5]);
 out(k,2:6) = NN / N;
end


for k=1:Nadd,
 for i=1:N,
     xi(i) = RV(temp(:,xi(i))');
 end
end
[NN XX] = histc(xi, [1 2 3 4 5]);

% Plotting and other aesthetics
figure(1); %clf;             % Create or clear Figure 1 if existing..
set(gcf, 'Color', [1 1 1]); % Change its default background color

subplot(3,2,1);             % Prepare to create 5 + 1 subplots
for k=1:5                   % For each of the first 5 subplots
 subplot(3,2,k);
 hold on;
 plot(out(:,1), out(:,1+k), 'r', 'LineWidth', 2);   % Plot P

 xlabel('t', 'FontSize', 25);           % X-label
 ylabel('{P(t)}', 'FontSize', 25);      % Y-label 
 
 L = legend(sprintf('%s', STATES{k}));  % Legend
 set(L, 'box', 'off');                  % No box around the legend
 
 set(gca, 'FontSize', 25, 'xlim', [0 dt*Nstep]);    % Change fonts
end

subplot(3,2,6);            % For the last subplot, represent the steady states
hold on;
B = bar(1:5, NN/N);
set(B, 'FaceColor', 'none', 'EdgeColor', 'r');

ylabel('{P_{\infty}}', 'FontSize', 25);
set(gca, 'XTick', 1:5, 'XTickLabel', STATES, 'xlim', [0 6], 'FontSize', 25);

subplot(3,2,1); ylim([0 1]);
subplot(3,2,2); ylim([0 1]);
subplot(3,2,3); ylim([0 1]);
subplot(3,2,4); ylim([0 0.06]);
subplot(3,2,5); ylim([0 1.5e-3]);
subplot(3,2,6); ylim([0 1]);
%set(gcf, 'units','normalized','outerposition',[0 0 1 1])
%print(gcf, '-dpng', '-opengl', 'analytical_solution.png');

disp('Steady-state occupancy probabilities:');
for k=1:5, disp(sprintf('P_%s = %.4f', STATES{k}, NN(k)/N)); end;
