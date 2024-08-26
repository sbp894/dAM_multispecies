clear;
clc;

fs= 20e3; % sampling frequency 
stimdur= 1;  % stimulus duration 
fc= 3e3; % The carrier 
dam_f_start= 2^4; % starting dAM frequency 
dam_f_end= 2^10.25; % end dAM frequency 

[stim, dam_traj_Hz, tStim] = create_dAM_stim(fs, stimdur, fc, dam_f_start, dam_f_end);

stim_env= abs(hilbert(stim)); 
stim_env= stim_env-mean(stim_env);


figure(4);
clf;

hold on;
subplot(211)
helper.plot_spectrogram(stim, fs, 40e-3, .95)
line(tStim*1e3, (fc+dam_traj_Hz)/1e3, 'color', 'r', 'linew', 2, 'linestyle', ':')
line(tStim*1e3, (fc+0*dam_traj_Hz)/1e3, 'color', 'r', 'linew', 2, 'linestyle', ':')
line(tStim*1e3, (fc-dam_traj_Hz)/1e3, 'color', 'r', 'linew', 2, 'linestyle', ':')
colorbar off;
ylim([(fc-dam_f_end*1.2) (fc+dam_f_end*1.2)]/1e3);
title(sprintf("Signal: F_c=%.1f kHz, F_m=%.0f to %.0f Hz", fc/1e3, dam_f_start, dam_f_end));

subplot(212)
hold on;
helper.plot_spectrogram(stim_env, fs, 40e-3, .95)
line(tStim*1e3, dam_traj_Hz/1e3, 'color', 'r', 'linew', 2, 'linestyle', ':')
colorbar off;
ylim([0 dam_f_end*1.2]/1e3);



