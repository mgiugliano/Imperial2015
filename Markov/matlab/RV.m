function out = RV(p)
 %
 % Simulate a discrete-valued random variable
 %
 % p: n x 1
 %
 % e.g., p = [0.2 0.4 0.4], returns values in {1 2 3}
 %
 % 2015 - Michele Giugliano, michele.giugliano@uantwerpen.be
 % http://www.uantwerpen.be/michele-giugliano

 out = -1;        % Output in case of an error
 
 d   = cumsum(p); % Now d = [0.2 0.6 1.0]
 r   = rand;      % Uniform random number in [0.0 1.0]
 
 h = 1;
 while h<=length(d),    % For each interval, let's check for it
  if (r <= d(h) ), out = h; return; end;
  h = h + 1;
 end
end