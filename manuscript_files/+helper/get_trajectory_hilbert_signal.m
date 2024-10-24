<<<<<<< HEAD
% function [fracSignal, filtSignal]= get_trajectory_signal(inSig, fs, freqTrajectory, filtParams)
function [fracSignal, filtSignal, sig_remod, sig_static_mod]= get_trajectory_hilbert_signal(inSig, fs, freqTrajectory, filtParams)

if nargin<3
    error('Need at least three inputs (input signal, sampling frequency, frequency trajectory along which we need to estimate power)');
end

if numel(inSig)~=numel(freqTrajectory)
    error('Length of inSig and freqTrajectory should be the same');
end

%% inputs should be column vectors
inSig= hilbert(inSig(:));
freqTrajectory= freqTrajectory(:);
freqTrajectory(isnan(freqTrajectory))= 0;
siglen= length(inSig);
stim_dur= siglen/fs;


if ~exist('filtParams', 'var') % doesn't exist
    filtParams= 2/stim_dur;
    d_lp = designfilt('lowpassiir','FilterOrder', 2, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
elseif isempty(filtParams) % exists but is empty
    filtParams= 2/stim_dur;
    d_lp = designfilt('lowpassiir','FilterOrder', 2, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
elseif isnumeric(filtParams) % exists and is numeric 
    d_lp = designfilt('lowpassiir','FilterOrder', 2, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
else % is a filter 
    d_lp= filtParams;
end

%%
phi_trajectory= -cumtrapz(freqTrajectory)/fs;

sig_demod_empirical= inSig .* exp(2*pi*1j*phi_trajectory);

filtSignal= filtfilt(d_lp, sig_demod_empirical);

sig_remod= real(filtSignal .* exp(-2*pi*1j*phi_trajectory));

sig_static_mod= real(filtSignal .* exp(-2*pi*1j*cumtrapz(mean(freqTrajectory))/fs));

filtSignal= abs(filtSignal); 
=======
% function [fracSignal, filtSignal]= get_trajectory_signal(inSig, fs, freqTrajectory, filtParams)
function [fracSignal, filtSignal, sig_remod, sig_static_mod]= get_trajectory_hilbert_signal(inSig, fs, freqTrajectory, filtParams)

if nargin<3
    error('Need at least three inputs (input signal, sampling frequency, frequency trajectory along which we need to estimate power)');
end

if numel(inSig)~=numel(freqTrajectory)
    error('Length of inSig and freqTrajectory should be the same');
end

%% inputs should be column vectors
inSig= hilbert(inSig(:));
freqTrajectory= freqTrajectory(:);
freqTrajectory(isnan(freqTrajectory))= 0;
siglen= length(inSig);
stim_dur= siglen/fs;


if ~exist('filtParams', 'var') % doesn't exist
    filtParams= 2/stim_dur;
    d_lp = designfilt('lowpassiir','FilterOrder', 2, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
elseif isempty(filtParams) % exists but is empty
    filtParams= 2/stim_dur;
    d_lp = designfilt('lowpassiir','FilterOrder', 2, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
elseif isnumeric(filtParams) % exists and is numeric 
    d_lp = designfilt('lowpassiir','FilterOrder', 2, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
else % is a filter 
    d_lp= filtParams;
end

%%
phi_trajectory= -cumtrapz(freqTrajectory)/fs;

sig_demod_empirical= inSig .* exp(2*pi*1j*phi_trajectory);

filtSignal= filtfilt(d_lp, sig_demod_empirical);

sig_remod= real(filtSignal .* exp(-2*pi*1j*phi_trajectory));

sig_static_mod= real(filtSignal .* exp(-2*pi*1j*cumtrapz(mean(freqTrajectory))/fs));

filtSignal= abs(filtSignal); 
>>>>>>> 033e1e8208ba0bc9e95b7a551303ab1999f126dc
fracSignal= filtSignal / rms(inSig); % normalize by total power