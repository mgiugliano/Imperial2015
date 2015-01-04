%
% exact_analytic_solution.m
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

Nstep = 500;    % Simulation time steps (exclusively used to plot traces)
Nadd  = 1000;   % Additional simulation time steps, until steady-state reached!

dt= 0.001;      % ms - iteration time step (it is an exact analytical method!)

a = 1;          % \alpha - [ms^-1] - values of the parameters
b = 7;          % \beta  - [ms^-1] - values of the parameters

STATES = {'A', 'B', 'C', 'D', 'E', 'F'};    % Aesthetics: name of states

% Initial condition ; note: sum(P) must be equal to 1 !!!!!!
%P = 0.2 * [1 1 1 1 1]';  
P = [0 0 1 0 0]';  

% Generator matrix of A <-> B <-> C <-> D <-> E, like "Kd" ion channels
Q = [-4*a     b         0      0     0; 
      4*a -(3*a+b)     2*b     0     0;
       0     3*a   -(2*a+2*b) 3*b    0;
       0      0        2*a -(a+3*b) 4*b;
       0      0         0      a   -4*b];
   
   
eQ = expm(Q * dt);  % Exponential of matrix Q * dt 
% see http://nl.mathworks.com/help/matlab/ref/expm.html
   

out = zeros(Nstep, 6);      % Output data structure, for plotting

for k=1:Nstep,		    % Main iteration loop, through time...
 % P = P + dt * Q * P;      % Forward Euler numerical method (NOT used)
 P = eQ * P;                % Iteration of the exact analytical solution

 out(k,1)   = k*dt;         % Log time data on the output structure
 out(k,2:6) = P';           % Log state-vars data on the output structure
end % end for

for k=1:Nadd,               % Iteration is continued until steady-state
 % P = P + dt * Q * P;      % Forward Euler numerical method (NOT used)
 P = eQ * P;                % Iteration of the exact analytical solution
end % end for


% Plotting and other aesthetics
figure(1); clf;             % Create or clear Figure 1 if existing..
set(gcf, 'Color', [1 1 1]); % Change its default background color

subplot(3,2,1);             % Prepare to create 5 + 1 subplots (i.e., P(t) for each state)
for k=1:5                   % For each of the first 5 subplots
 subplot(3,2,k);
 plot(out(:,1), out(:,1+k), 'k', 'LineWidth', 2);   % Plot P(t), for each of the 5 states

 xlabel('t',      'FontSize', 25);      % X-label
 ylabel('{P(t)}', 'FontSize', 25);      % Y-label 
 
 L = legend(sprintf('%s', STATES{k}));  % Legend, indicating the name of the corresponding states
 set(L, 'box', 'off');                  % No box around the legend
 
 set(gca, 'FontSize', 25, 'xlim', [0 dt*Nstep]);    % Change fonts for the axes ticks
end

subplot(3,2,6);            % For the last subplot, let's represent the steady states values
B = bar(1:5, P);	   % Plot the current value of P as a bar plot (i.e., at time t = t_end)
set(B, 'FaceColor', 'k', 'EdgeColor', 'k');	% Make the bar plot nicely colored

ylabel('{P_{\infty}}', 'FontSize', 25);		% Y-label
set(gca, 'XTick', 1:5, 'XTickLabel', STATES, 'xlim', [0 6], 'FontSize', 25);

%set(gcf, 'units','normalized','outerposition',[0 0 1 1])	% Unused, does not work nicely
%print(gcf, '-dpng', '-opengl', 'analytical_solution.png');	% Unused, does not work nicely

disp('Steady-state occupancy probabilities:');			% Print the steady-state values
for k=1:5, disp(sprintf('P_%s = %.4f', STATES{k}, P(k))); end;	% on the command window
