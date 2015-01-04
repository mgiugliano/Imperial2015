function [T varargout] = readfiletocells(filename)
% READFILETOCELLS reads a text file in which every line contains the times
% at which a neuron fired.
%
% T = readfiletocells(filename)
% [T N] = readfiletocells(filename)
% [T N Tlim] = readfiletocells(filename)
%
% Arguments:
%   filename - the name of the file.
% 
% Returns:
%   T - an array of cells containing the spike times. The number of cells
%   is equal to the number of lines in the file.
%   N - an array containing the number of spikes fired by each neuron.
%   Tlim - an array containing the times of the last spike fired by each
%   neuron.
% 

%
%   Author: Daniele Linaro - August 2009
%

fid = fopen(filename,'r');
T = {};
row = 1;
l = fgetl(fid);
while l ~= (-1)
    [t,r] = strtok(l);
    if isempty(t)
        break;
    end
    ind = 1;
    T{row}(ind) = str2double(t);
    ind = ind + 1;
    while ~ isempty(r)
        [t,r] = strtok(r);
        if isempty(t)
            break;
        end
        T{row}(ind) = str2double(t);
        ind = ind + 1;
    end
    row = row+1;
    l = fgetl(fid);
end
fclose(fid);
T = T(:);

if nargout == 2
    N = zeros(length(T),1);
    for ii=1:length(T)
        N(ii) = length(T{ii});
    end
    varargout{1} = N;
end

if nargout == 3
    Tlim = [1e20 -1e20];
    for ii=1:length(T)
        if T{ii}(1) < Tlim(1)
            Tlim(1) = T{ii}(1);
        end
        if T{ii}(end) > Tlim(2)
            Tlim(2) = T{ii}(end);
        end
    end
    varargout{2} = Tlim;
end

