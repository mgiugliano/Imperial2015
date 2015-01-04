function rasterplot(spktimes,varargin)
% RASTERPLOT plots a raster plot of spike times.
% 
% rasterplot(spktimes,...)
% 
% Arguments:
%   spktimes - the times at which spikes were fired. It must be an array
%   of cells, where each cell contains the spike times of a single neuron
%   or trial.
% 
% Additional (optional) arguments will be passed unchanged to the function 
% plot.
% 

% 
% Author: Daniele Linaro - September 2009.
% 

hold on;
for ii=1:length(spktimes)
    for jj=1:length(spktimes{ii})
        plot([spktimes{ii}(jj),spktimes{ii}(jj)],[ii-1,ii],varargin{:});
    end
end
