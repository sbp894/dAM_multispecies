clear;
clc;

% Two key functions: both are used in this script (with simulated FFRs) 
% 1) get_trajectory_hilbert_signal(signal, sampling_frequency, desired_trajectory, half_bandwidth_for_filtering)
% 2) get_freq_trajectory_power(signal, sampling_frequency, desired_trajectory, flag_plotPSD, half_bandwidth_for_filtering)
%   if flag_plotPSD = true, it plots the power spectral density (PSD)

%% Create a stimulus with a desired trajectory 
fs_stim= 12e3; % sampling frequency 
stimdur= 1;  % stimulus duration 
stimlen= round(fs_stim*stimdur); % stimulus length in samples 
dAM_f_start= 2^4; % starting dAM frequency 
dAM_f_end= 2^10.25; % end dAM frequency 

fc= 3e3; % The carrier 
mod_depth= 1; % modulation depth 
t_stim= (1:stimlen)/stimlen; % time 

dAM_stim_traj_Hz= 2.^linspace(log2(dAM_f_start), log2(dAM_f_end), stimlen); % Desired dAM trajectory in stimulus 

phi_instantaneous= -cumtrapz(dAM_stim_traj_Hz)/fs_stim; % Formula for phase (in terms of frequency)
stimMod= cos(2*pi*phi_instantaneous); % modulator 
stimCar= sin(2*pi*fc*t_stim); % carrier 
dAM_stim= (1-mod_depth*stimMod(:))/2.*stimCar(:); % dAM stimulus 

%% Read (or simulate) data
fs_resp= 5e3;
ffr_env_resp= abs(hilbert(abs(dAM_stim))); % treat the stimulus envelope as the FFR response, rectifying the stimulus to introduce harmonics 
ffr_env_resp= ffr_env_resp + .25*randn(size(ffr_env_resp)); % add a bit of noise 
ffr_env_resp= ffr_env_resp-mean(ffr_env_resp);  % subtract the mean 
ffr_env_resp= resample(ffr_env_resp, fs_resp, fs_stim); % downsample because generally FFRs have a lower frequency than stimuli 
ffr_env_resp= ffr_env_resp/max(abs(ffr_env_resp));
t_resp= (1:length(ffr_env_resp))/fs_resp;

dAM_resp_traj_Hz= interp1(t_stim(:), dAM_stim_traj_Hz(:), t_resp(:));  % Desired dAM trajectory in the response (which has a different fs) 

%% Get time-domain filtered output along the desired trajectory 
Filter_HalfWindow_Hz= 12;
[dAM_pow_frac_tracked, dAM_pow_total_tracked]= get_trajectory_hilbert_signal(ffr_env_resp, fs_resp, dAM_resp_traj_Hz, Filter_HalfWindow_Hz);

%% Get frequency-domain power along the desired trajectory 
plotPSD = 0;
[fractional_Power, total_Power, PSD, Freq_psd]= get_freq_trajectory_power(ffr_env_resp, fs_resp, dAM_resp_traj_Hz, plotPSD, Filter_HalfWindow_Hz);

%% plot everything 
figure(1);
clf;
set(gcf, "Units", "inches", 'Position', [1, 1, 10, 3])

% plot the stimulus 
ax(1)= subplot(3, 5, 1:3);
plot_dAM_spectrogram(dAM_stim, fs_stim);
hold on
box off;
plot(t_stim*1e3, 1+max(ylim())+dAM_stim)
ylim([0, 8])
xlabel('')
ylabel('Stimulus')


% plot the response 
ax(2)= subplot(3, 5, 6:8);
plot_dAM_spectrogram(ffr_env_resp, fs_resp);
hold on
box off;
plot(t_resp*1e3, .5+max(ylim())+.5*ffr_env_resp)
ylim([0, 3.5])
xlabel('')
ylabel('FFR')

% show filtered signal along the desired trajectory 
ax(3)= subplot(3, 5, 11:13);
hold on;
plot(t_resp, dAM_pow_frac_tracked, 'DisplayName', "Fractional");
plot(t_resp, dAM_pow_total_tracked, 'DisplayName', "Absolute");
legend('show', 'box', 'off')
ylim([0, 1.1])
xlablhan3= xlabel('Time (s)', 'Units', 'normalized');
xlablhan3.Position(2)= -.2;
ylabel('Filtered sig')

% plot demodulated spectrum (transposed to a random frequency value, 140 Hz in this
% case) to show how dAM is mapped onto a single frequency 
ax(4)= subplot(2, 5, [4, 5, 9, 10]);
plot(Freq_psd/1e3, PSD);
xlabhan= xlabel('Freq (kHz)', 'Units','normalized');
xlabhan.Position(2)= -.06;
ylabel('Power (dB)')
set(gca, 'XScale', 'log')
xlim([50, fs_resp/2]/1e3);
ylim([-60, -10]);
box off;
title('demodulated spectrum');

set(gca, 'XTick', [.1, .2, .5, 1, 2, 4, 8])

set(findall(gcf, "-property", "FontSize"), "FontSize", 9);

%% 
xc= .04; 
xs= .055; 
xw1= .6; 
xw2= 1-xs-xc-xw1-.01;
yc= .11;
ys= .07;
yw1= .245;
yw2= yc + 2*ys+3*yw1 -.15;

set(ax(3), 'Units', 'normalized', 'Position', [xc, yc, xw1, yw1]);
set(ax(2), 'Units', 'normalized', 'Position', [xc, yc+yw1+ys, xw1, yw1]);
set(ax(1), 'Units', 'normalized', 'Position', [xc, yc+2*yw1+2*ys, xw1, yw1]);

set(ax(4), 'Units', 'normalized', 'Position', [xc+xw1+xs, yc, xw2, yw2]);
