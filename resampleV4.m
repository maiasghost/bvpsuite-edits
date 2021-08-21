function tosignal = resampleV4(fromgrid,fromsignal,togrid)
% Oversampling that guarantees that tosignal remains positive
% Written by GS, TU Wien, 24 August 2006
% V4 is original version

signal = abs(fromsignal);                   % For robustness only
meansig = norm(signal,1)/length(signal);
mag = max(signal);
signal = signal/mag;                        % Normalize to maximum of 1
signal = signal + 1e-3*exp(-5*signal);      % Limiter lifts values near zero 1e-3
signal = log(signal);                       % Transform to logarithmic scale
hi = max(signal);
lo = min(signal);
newsignal = spline(fromgrid,signal,togrid); % Resample on new grid
signal = max(lo,min(hi,newsignal));         % Make sure amplitude doesn't grow
signal = exp(newsignal);                    % Now mapped to new grid and positive
signal = signal + 1e-5*exp(-3*signal);      % Limiter lifts values near zero 1e-2
signal = bcTDFlogV4(signal);                % Smoothing
signal = mag*signal/max(signal);            % Restore signal amplitude
newmean = norm(signal,1)/length(signal);    % Correct the signal's mean value
signal = signal + meansig - newmean;
signal = max(0,signal);                     % Protect against negative values
signal = signal + 1e-5*exp(-3*signal);      % Limiter lifts values near zero 1e-3
tosignal = signal;                          % Output