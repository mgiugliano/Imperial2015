%
% variable_step.m
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
N = 5000;       % Number of independent and identical channels

Nstep = 500000;    % Simulation time steps (exclusively used to plot traces)
Nadd  = 1000000;   % Additional simulation time steps, until steady-state reached!

dt= 0.000001;      % ms - simulation time step 

a = 1;          % \alpha - [ms^-1] - values of the parameters
b = 7;          % \beta  - [ms^-1] - values of the parameters

STATES = {'A', 'B', 'C', 'D', 'E', 'F'};    % Aesthetics: name of states

% Equivalent initial condition 
%xi = ceil(rand(N,1)*5); % equivalent to P = [0 0 1 0 0]';, xi \in {1,2,3,4,5} 
xi = 3 * ones(N,1); % equivalent to P = [0 0 1 0 0]';, xi \in {1,2,3,4,5} 

% Generator matrix of A <-> B <-> C <-> D <-> E, like "Kd" ion channels
Q = [-4*a     b         0      0     0; 
      4*a -(3*a+b)     2*b     0     0;
       0     3*a   -(2*a+2*b) 3*b    0;
       0      0        2*a -(a+3*b) 4*b;
       0      0         0      a   -4*b];
      
temp = eye(5) + Q*dt;        % Useful temporary matrix
d   = -diag(Q);              % Diagonal of Q
tmq = Q .* (ones(5) - eye(5));
for h=1:5,
    tmq(:,h) = tmq(:,h) / sum(tmq(:,h));
end

out = zeros(Nstep, 6);      % Output data structure, for plotting 

k = 1;
while k<=Nstep,
 T   = -(1./d(xi)) .* log(rand(N,1));       % Random state life-times
 mT  = min(T);                              % The minimum is the one we want to select..
 idx = find(T == mT);                       % Find the specific channel who made the transition..
 
 xi(idx) = RV(tmq(:,xi(idx))');             % Now, given that a transition occurred, which one is it?
 
 ndt = ceil(mT/dt);                         % How much time passed since the last transition?
 
 L = min([k+ndt, Nstep]);
 out(k:L,1)   = (k:L)*dt;
 [NN XX] = histc(xi, [1 2 3 4 5]);
 for h=1:5,
     out(k:L,h+1) = NN(h) / N;
 end
 k = L + 1;
end


k = 1;
while k<=Nadd,
    
 T   = -(1./d(xi)) .* log(rand(N,1));       % Random state life-times
 mT  = min(T);                              % The minimum is the one we want to select..
 idx = find(T == mT);                       % Find the specific channel who made the transition..
 
 xi(idx) = RV(tmq(:,xi(idx))');             % Now, given that a transition occurred, which one is it?
 
 ndt = ceil(mT/dt);                         % How much time passed since the last transition?
 L = min([k+ndt, Nadd]);
 k = L + 1;
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
