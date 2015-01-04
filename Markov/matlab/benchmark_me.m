%
% benchmark_me.m
% (you must slightly alter both variable_step.m and uniform_step.m, as indicated below)
%
% (c) 2015 - Michele Giugliano, michele.giugliano@uantwerpen.be
% http://www.uantwerpen.be/michele-giugliano
%
% Reference: Linaro D, Giugliano M (2014) Markov models of ion channels.
% Springer Encyclopedia of Computational Neuroscience
% doi:10.1007/978-1-4614-7320-6_131-1
%

clear all;      % Clear all varibles from the memory
close all;      % Close all files and figures
clc;            % Clear the command window

disp('Important: before running this script, ensure that, in variable_step.m and uniform_step.m, the line corresponding to');
disp('N = 100;       % Number of independent and identical channels');
disp('is commented!');
return;		% Remove or comment this line to launch the benchmarking..

Nrange   = [1 2 5 8 10 20 50 80 100 200 500 800 1000 2000];
Nrep     = 10;

%% CPU time probing for 'uniform_step.m'
disp('Warning: this might take 1h or more to complete!');
tt = zeros(length(Nrange), 3);		% Data structure to store CPU times, for increasing N
for kkk=1:length(Nrange),		% Main loop over the values of N (Nrange)
    N = Nrange(kkk);       		% Number of independent and identical channels
    tmp_toc = zeros(Nrep,1);		% Temporary structure, initialised here
    for r=1:Nrep,			% The CPU execution time is sampled for Nrep times, not once!
        close all;			% All figures and files are closed
        tic;				% The stopwatch is started
        %variable_step;
        uniform_step;
        tmp_toc(r) = toc;		% The stopwatch is stopped
    end
    tt(kkk,1) = N			% The output data structure is filled with the relevant information
    tt(kkk,2) = mean(tmp_toc);		% such as the average CPU time and
    tt(kkk,3) = std(tmp_toc);		% the standard deviation (currently not represented graphically)
end
save('benchmark_uniform_step.mat', 'tt'); % These data are stored on disk!

%% CPU time probing for 'variable_step.m'
tt = zeros(length(Nrange), 3);          % Data structure to store CPU times, for increasing N
for kkk=1:length(Nrange),               % Main loop over the values of N (Nrange)
    N = Nrange(kkk);                    % Number of independent and identical channels
    tmp_toc = zeros(Nrep,1);            % Temporary structure, initialised here
    for r=1:Nrep,                       % The CPU execution time is sampled for Nrep times, not once!
        close all;                      % All figures and files are closed
        tic;                            % The stopwatch is started
        variable_step;
        %uniform_step;
        tmp_toc(r) = toc;               % The stopwatch is stopped
    end
    tt(kkk,1) = N                       % The output data structure is filled with the relevant information
    tt(kkk,2) = mean(tmp_toc);          % such as the average CPU time and
    tt(kkk,3) = std(tmp_toc);           % the standard deviation (currently not represented graphically)
end
save('benchmark_variable_step.mat', 'tt');  % These data are stored on disk!


%% Let's plot the benchmark results
clear all;      % Clear all varibles from the memory
close all;      % Close all files and figures
clc;            % Clear the command window

uniform = load('benchmark_uniform_step.mat');	% Results are loaded from disk...
variable= load('benchmark_variable_step.mat');	% Results are loaded from disk...

figure(2); clf;             % Create or clear Figure 1 if existing..
set(gcf, 'Color', [1 1 1]); % Change its default background color

P = plot(uniform.tt(:,1), uniform.tt(:,2), 'ko-', variable.tt(:,1), variable.tt(:,2), 'rp-');
set(P, 'LineWidth', 2, 'MarkerSize', 10);% Aesthetics..
ylim([-10 300])				 % Adjust the vertical limits to enhance differences in CPU times
xlabel('N', 	       'FontSize', 25);  % X-label
ylabel('CPU time [s]', 'FontSize', 25);	 % Y-label
set(gca, 'FontSize', 20);  		 % Aesthetics..
